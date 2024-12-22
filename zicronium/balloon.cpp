#define BALLOON_DEBUG

/**
 *
 * BALLOON - C++ lightweight zlib implementation
 *
 * Program written by muffinshades
 *
 * Copyright (c) 2024 muffinshades
 *
 * You can do what ever you want with the software but you must
 * credit the author (muffinshades) with the original creation
 * of the software. Idk what else to put here lmao.
 *
 * Balloon Notes:
 *
 * This library is a implementation of the zlib or Inflate / Deflate
 * compression algorithm.
 *
 * Right now compression speeds are around 4-6mb/s and decompression
 * speeds are much faster. This isn't the worlds fastest implementation,
 * but its decently fast and lightweight. One day I will improve the LZ77
 * hash functions, but for now it's gonna stay at around 5mb/s.
 *
 * This program should be able to function with any other inflate / deflate
 * implementation apart from the whole compression level calculations being
 * different since I didn't entirley implement lazy and good matches into
 * the lz77 functions. I also didn't add a whole fast version for everything
 * since this is a relativley light weight library. One day I do plan on adding
 * these functions and making a even better implementation of zlib.
 *
 */

#include "balloon.hpp"
#include "ByteStream.hpp"
#include "linked_map.hpp"
#include "micro_vector.h"
#include <algorithm>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <queue>


 //bin
struct bin {
    byte* dat;
    size_t sz;
};

//MACROS and constants
#define INRANGE(v,m,x) ((v) >= (m) && (v) < (x))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))
#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define WINDOW_BITS 0xf
#define MIN_MATCH 3
#define MEM_LEVEL 8

#ifndef INT_MAX
#define INT_MAX 0xffffffff
#endif

#define LZ77_MIN_MATCH 0x003 // min match size for back reference
#define LZ77_MAX_MATCH 0x102

//from zlib on github
const int compression_table_old[10][5] = {
    /*      good lazy nice max_chain_length */
    /* 0 */ {0,    0,  0,    0},    /* store only */
    /* 1 */ {4,    4,  8,    4},    /* max speed, no lazy matches */
    /* 2 */ {4,    5, 16,    8},
    /* 3 */ {4,    6, 32,   32},
    /* 4 */ {4,    4, 16,   16},    /* lazy matches */
    /* 5 */ {8,   16, 32,   32},
    /* 6 */ {8,   16, 128, 128},
    /* 7 */ {8,   32, 128, 256},
    /* 8 */ {32, 128, 258, 1024},
    /* 9 */ {32, 258, 258, 4096}    /* max compression */
};

//lz77 constants
const i32 LengthExtraBits[] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0 };
const i32 LengthBase[] = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258 };
const i32 DistanceExtraBits[] = { 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13 };
const i32 DistanceBase[] = { 1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577 };

const uint32_t MOD_ADLER = 65521; //0xfff1
constexpr u32 ADLER_MASK = 0xffff;

#define ADLER32_BASE 1

void adler32_compute_next(u32& cur, const byte val) {
    u32 a = cur & ADLER_MASK,
        b = (cur >> 0x10) & ADLER_MASK;

    a = (a + val) % MOD_ADLER;
    b = (b + a) % MOD_ADLER;

    cur = (b << 0x10) | a;
}

#define RLE_Z2_BASE 0xb
#define RLE_Z1_BASE 0x3

#define RLE_Z2_MASK 0x7f
#define RLE_Z1_MASK 0x7

#define RLE_L_BITS  0x2
#define RLE_Z1_BITS 0x3
#define RLE_Z2_BITS 0x7

#define RLE_L_MASK 0x3
#define RLE_L_BASE 0x3

/**
 *
 * BitStream Functions
 *
 * Deflate code starts after this
 *
 */

BitStream::BitStream(u32* bytes, size_t len) {
    this->bytes = new u32[len];
    memcpy(this->bytes, bytes, len * sizeof(u32));
    this->bsz = len * 8;
    this->sz = len;
    this->asz = sz;
    this->rPos = this->sz;
    this->pos = 0;
}

BitStream::BitStream(size_t len) {
    assert(len > 0);
    this->bytes = new u32[len];
    memset(this->bytes, 0, len * sizeof(u32));
    this->bsz = len * 8;
    this->asz = len;
    this->rPos = 0;
    this->sz = 0;
}

i32 BitStream::readBit() {
    if (this->lBit <= 0) {
        assert(this->pos < this->sz && this->pos >= 0);

        //advance a byte
        this->cByte = this->bytes[this->pos];
        this->pos++;
        this->lBit = 8;
    }

    this->lBit--;

    //extract le bit
    i32 bit = this->cByte & 1; //extract first bit
    this->cByte >>= 1; //advance a bit
    return bit; //return extracted bit
}

i32 BitStream::readByte() {
    this->lBit = 0;
    return this->bytes[this->pos++];
}

i32 BitStream::readNBits(i32 nBits) {
    u32 r = 0, c = 0;

    while (c < nBits) {
        i32 bit = this->readBit() & 1;
        r |= (bit << (c++));
    }
    return r;
}

i32 BitStream::tell() {
    return this->pos;
}

i32 BitStream::tellw() {
    return this->rPos;
}

void BitStream::seekw(i32 wPos) {
    assert(wPos >= 0 && wPos < this->asz);
    this->rPos = wPos;
}

void BitStream::seek(i32 pos) {
    assert(pos >= 0 && pos < this->asz);
    this->pos = pos;
}

void BitStream::checkWritePosition() {
    if (this->rPos >= this->asz)
        this->allocNewChunk();
}

i32 BitStream::readValue(size_t size) {
    i32 r = 0;
    for (size_t i = 0; i < size; i++)
        r |= (this->readByte() << (i * 8));

    return r;
}

u32* BitStream::bitsToBytes(u32* bits, size_t len) {
    u32* res = new u32[(u32)ceil((float)len / 8.0f)];

    i32 bCollect = 0, iBit = 0, cByte = 0;
    for (i32 i = len - 1; i >= 0; i--) {
        bCollect |= (bits[i]) << iBit;
        iBit++;
        if (iBit >= 8) {
            res[cByte++] = bCollect;
            bCollect = 0;
            iBit = 0;
        }
    }

    if (iBit > 0) res[cByte] = bCollect;

    return res;
}

u32* BitStream::bitsToBytes(u32 bits, size_t len) {
    u32* res = new u32[(u32)ceil((float)len / 8.0f)];

    i32 bCollect = 0, iBit = 0, cByte = 0;
    for (i32 i = len - 1; i >= 0; i--) {
        bCollect |= ((bits >> i) & 1) << iBit;
        iBit++;
        if (iBit >= 8) {
            res[cByte++] = bCollect;
            bCollect = 0;
            iBit = 0;
        }
    }

    if (iBit > 0) res[cByte] = bCollect;

    return res;
}

void BitStream::writeBit(u32 bit) {
    //now write bit
    this->bytes[this->rPos] |= ((bit & 1) << (this->rBit++)); //& with 1 to just get first bit
    //advance a byte if were out of range
    if (this->rBit >= 8) {
        this->rBit = 0;
        this->rPos++;
        this->sz++;
        this->checkWritePosition();
    }
}

void BitStream::writeNBits(u32* bits, size_t len) {
    for (size_t i = 0; i < len; i++)
        this->writeBit(bits[i]);
}

template<typename _T> void t_writeValue(BitStream* s, _T val) {
    size_t vsz = sizeof(_T) * 8;

    for (i32 i = 0; i < vsz; i++)
        s->writeBit((val >> i) & 1);
}

//dll template fix
void BitStream::writeValue(byte val) {
    t_writeValue(this, val);
}

void BitStream::writeValue(short val) {
    t_writeValue(this, val);
}

void BitStream::writeValue(i32 val) {
    t_writeValue(this, val);
}

void BitStream::writeValue(u32 val) {
    t_writeValue(this, val);
}

void BitStream::writeValue(i64 val) {
    t_writeValue(this, val);
}

void BitStream::writeValue(u64 val) {
    t_writeValue(this, val);
}

//
void BitStream::allocNewChunk() {
    this->asz += 0xffff;
    u32* tBytes = new u32[this->asz];
    memset(tBytes, 0, sizeof(u32) * this->asz);
    memcpy(tBytes, this->bytes, sizeof(u32) * this->sz);
    delete[] this->bytes;
    this->bytes = tBytes;
}

void BitStream::clip() {
    //wtf???
    u32* tBytes = new u32[this->asz];
    u32 osz = this->asz;
    memcpy(tBytes, this->bytes, this->asz);
    this->asz = this->sz;
    delete[] this->bytes;
    this->bytes = new u32[this->asz];
    memset(this->bytes, 0, sizeof(u32) * this->asz);
    memcpy(this->bytes, tBytes, osz);
    delete[] tBytes;
}

void BitStream::calloc(size_t sz) {
    if (this->bytes)
        delete[] this->bytes;

    this->bytes = new u32[sz];

    this->asz = sz;
    this->bsz = sz * 8;
    this->pos = this->rPos = this->lBit = this->cByte = this->sz = 0;

    i32 pos = 0, rPos = 0;
    i32 lBit = 0;
    u32 cByte = 0;
    u32 rBit = 8;
}


struct huffman_node {
    huffman_node* left = nullptr, * right = nullptr; //left and right pointers
    u32 val = 0;
    size_t count = 0, depth = 0, cmpCount = 0;
    bool leaf = false;
};

//
struct treeComparison {
    bool operator()(huffman_node* a, huffman_node* b) {
        if (a->cmpCount == b->cmpCount)
            return a->depth > b->depth;
        return a->cmpCount > b->cmpCount;
    }
};

/**
 *
 * GenTreeFromCounts
 *
 * Generate a huffman tree from a bunch of char counts
 * Will auto filter out any characters with a count of 0
 *
 * Params:
 *  template: [size_t] alphaSz - size of the alphat (number of characters in count array)
 *  [size_t*] counts[alphaSz] - array containing the character counts in the given message
 *
 * Return:
 *  huffman_node * - root of the constructed huffman tree
 *
 * Errors:
 *  will return nullptr if counts are in valid in anyway
 *
 */

template<size_t alphaSz> huffman_node* GenTreeFromCounts(size_t* counts, i32 maxLen = -1) {
    if (counts == nullptr || alphaSz <= 0)
        return nullptr;

    //convert all of the char counts into a priority queue of huffman nodes
    std::priority_queue<huffman_node*, std::vector<huffman_node*>, treeComparison> tNodes;
    size_t _val = 0;
    for (size_t* count = counts, *end = counts + sizeof(u32) * alphaSz;
        end - count > 0;
        count += sizeof(u32)
        )
        if (*count > 0)
            tNodes.push(new huffman_node{
                .val = (u32) _val++,
                .count = *count,
                .cmpCount = *count,
                .leaf = true
                });


    //tree root
    huffman_node* root = nullptr;

    //construct tree
    while (tNodes.size() > 1) {
        huffman_node* newNode = new huffman_node;
        newNode->left = tNodes.top();
        tNodes.pop();
        newNode->right = tNodes.top();
        tNodes.pop();
        newNode->depth = MAX(newNode->left->depth, newNode->right->depth) + 1;
        newNode->cmpCount = newNode->count = newNode->left->count + newNode->right->count;

        //max length thingy for le tree
        if (maxLen >= 0 && newNode->depth >= maxLen - 1) newNode->cmpCount = INT_MAX;

        root = newNode;
        tNodes.push(newNode);
    }

    return root;
}

/**
 * 
 * getTreeBitLens
 * 
 * gets the length of the bit lengths for a tree
 * 
 */
size_t* getTreeBitLens(huffman_node* tree, size_t alphabetSize, i32 currentLen = 0, size_t* cbl = nullptr) {
    if (cbl == nullptr) {
        cbl = new size_t[alphabetSize];
        ZeroMem(cbl, alphabetSize);
    }

    if (tree->right || tree->left) {
        if (tree->left) getTreeBitLens(tree->left, alphabetSize, currentLen + 1, cbl);
        if (tree->right) getTreeBitLens(tree->right, alphabetSize, currentLen + 1, cbl);
    }
    else {
        if (tree->val < alphabetSize)
            //cbl[tree->val] = currentLen > 0 ? currentLen : 1;
            cbl[tree->val] = max(currentLen, 1);
        else
            std::cout << "Invalid bl_symbol [" << tree->val << "]" << "  " << alphabetSize << std::endl;
    }

    return cbl;
}

/**
 *  
 * getBitLenCounts
 * 
 * gets the bit lengths counts for a tree uh yeah
 * 
 */
size_t* getBitLenCounts(size_t* bitLens, size_t blLen, size_t MAX_BITS) {
    size_t* res = new size_t[MAX_BITS + 1];
    ZeroMem(res, MAX_BITS + 1);

    for (i32 i = 0; i < blLen; i++)
        res[bitLens[i]]++;

    return res;
}

/**
 * 
 * BitLengthsToHTree
 * 
 * Converts the bit lengths of a tree into a canonical tree
 * 
 * bitLens - bitLengths
 * blLEn - number of bit lengths
 * 
 */
huffman_node* BitLengthsToHTree(size_t* bitLens, size_t blLen, size_t aLen) {
    const size_t MAX_BITS = ArrMax<size_t>(bitLens, blLen);
    size_t* blCount = getBitLenCounts(bitLens, blLen, MAX_BITS);

    std::vector<i32> nextCode = { 0,0 };

    //get the next code for the tree items
    for (i32 codeIdx = 2; codeIdx <= MAX_BITS; codeIdx++)
        nextCode.push_back((nextCode[codeIdx - 1] + blCount[codeIdx - 1]) << 1);

    //construct tree
    i32 i = 0;
    huffman_node* res = new huffman_node();

    while (i < aLen && i < blLen) {
        if (bitLens[i] != 0) {
            InsertNode(nextCode[bitLens[i]], bitLens[i], i, res);
            nextCode[bitLens[i]]++;
        }
        i++;
    }

    //free and return
    delete[] blCount;

    return res;
}

/**
 *
 *  CovnertTreeToCanonical
 * 
 * converts a tree into its canonical form :O
 * 
 */
huffman_node* CovnertTreeToCanonical(huffman_node* tree, size_t alphabetSize, bool deleteOld = false) {
    //get bit lengths of each code
    size_t* bitLens = getTreeBitLens(tree, alphabetSize);

    //start to create le tree
    huffman_node* cTree = BitLengthsToHTree(bitLens, alphabetSize, alphabetSize);

    //free memory and return new tree
    delete[] bitLens;

    if (deleteOld)
        TreeFree(tree);

    return cTree;
}

//helper function that increase bit length for all child nodes of a node
void _bl_inc(size_t** bl, huffman_node* node, const size_t charMax) {
    if (!bl || !node) return;

    if (node->leaf)
        if (node->val < charMax)
            *bl[node->val]++; //increase bit length
    else {
        if (node->left)
            _bl_inc(bl, node->left, charMax);
        if (node->right)
            _bl_inc(bl, node->right, charMax);
    }
}

/**
 * 
 * GenCanonicalTreeFromCounts
 * 
 * generates a canonical tree slight faster then 
 * creating a temp tree and them using convert
 * tree to canonical.
 * 
 * Basically all it does is create a temp tree but
 * whilst creating it, bit lengths are counted so
 * you don't need to do the slow process of
 * also counting bit lengths
 * 
 */
template<size_t alphaSz> huffman_node* GenCanonicalTreeFromCounts(size_t* counts, i32 maxLen = -1) {
    if (counts == nullptr || alphaSz <= 0)
        return nullptr;

    //convert all of the char counts into a priority queue of huffman nodes
    std::priority_queue<huffman_node*, std::vector<huffman_node*>, treeComparison> tNodes;
    size_t _val = 0;
    for (size_t* count = counts, *end = counts + sizeof(u32) * alphaSz;
        end - count > 0;
        count += sizeof(u32)
    )
        if (*count > 0)
            tNodes.push(new huffman_node{
                .val = (u32)_val++,
                .count = *count,
                .cmpCount = *count,
                .leaf = true
            });


    //tree root
    huffman_node* root = nullptr;

    size_t* bitLens = new size_t[alphaSz];
    ZeroMem(bitLens);

    //construct tree
    while (tNodes.size() > 1) {
        huffman_node* newNode = new huffman_node;
        newNode->left = tNodes.top();
        tNodes.pop();
        newNode->right = tNodes.top();
        tNodes.pop();
        newNode->depth = MAX(newNode->left->depth, newNode->right->depth) + 1;
        newNode->cmpCount = newNode->count = newNode->left->count + newNode->right->count;

        //max length thingy for le tree
        if (maxLen >= 0 && newNode->depth >= maxLen - 1) newNode->cmpCount = INT_MAX;

        root = newNode;
        tNodes.push(newNode);

        //as we construct tree keep track of bitlengths
        _bl_inc(&bitLens, newNode, alphaSz);
    }

    TreeFree(root); //delete temporary tree

    //generate canonical tree
    return BitLengthsToHTree(bitLens, alphaSz, alphaSz);
}

/**
 *
 * GenCharacterCounts
 *
 * gets the count of characters in a bunch of bytes
 *
 * Params:
 *  [byte*] data - data to get the char counts of
 *  [size_t] sz - size of the data in bytes
 *
 * Return:
 *  u32* - char counts, will have a size of 0xff or BYTE_MAX
 *
 */

u32* GetCharacterCounts(byte* data, size_t sz) {
    u32* rCounts = new u32[0xff];
    ZeroMem(rCounts, 0xff);
    byte* d_ptr = data;
    while (sz--) rCounts[*d_ptr++]++;
    return rCounts;
};

//same thing as prev function just for u32* instead of byte*
u32* GetCharacterCounts(u32* data, size_t sz) {
    u32* rCounts = new u32[0xff];
    ZeroMem(rCounts, 0xff);
    u32* d_ptr = data;
    while (sz--) rCounts[*d_ptr++]++;
    return rCounts;
};


/**
 *
 * InsertNode
 *
 * Function used to insert values into a huffman tree during
 * the Inflate process. This is to construct the tree from a
 * list of codes generated by blListToTree
 *
 */

void InsertNode(u32 code, size_t codeLen, u32 val, huffman_node* tree) {
    huffman_node* current = tree;

    for (i32 i = codeLen - 1; i >= 0; i--) {
        u32 bit = (code >> i) & 1;

        if (!!bit) {
            if (current->right == nullptr)
                current->right = new huffman_node{
                    .count = 0
            };

            current = current->right;
        }
        else {
            if (current->left == nullptr)
                current->left = new huffman_node{
                    .count = 0
            };

            current = current->left;
        }
    }

    current->val = val;
}

/**
 *
 * TreeFree
 *
 * Used to free huffman trees
 *
 */
void TreeFree(huffman_node* root) {
    if (root == nullptr) return;

    if (root->left) TreeFree(root->left);
    if (root->right) TreeFree(root->right);

    delete root;
    root = nullptr;
}

//some lz77 helper functions
size_t lz77_get_len_idx(size_t len) {
    size_t lastIdx = 0;

    for (i32 l : LengthBase) {
        if (l > len) break;

        lastIdx++;
    }

    return lastIdx;
}

//get codes and indexs from distances
i32 lz77_get_dist_idx(size_t dist) {
    size_t lastIdx = 0;

    for (i32 l : DistanceBase) {
        if (l > dist) break;

        lastIdx++;
    }

    return lastIdx;
}

/**
 *
 * lz77 encode
 *
 * function for the lz77 compression algorithm
 *
 */

//window ref
struct _lz_win_ref {
    size_t len = 0;
    byte* ptr = nullptr;
};

/**
 * 
 * lzr
 * 
 * result from lz77 encoding
 * 
 */
struct lzr {
    u32* encDat;
    size_t encSz, winBits;
    huffman_node* distTree;
    u32* charCounts;
};

/**
 *
 * lz_inst
 * 
 * internal instance thing used in lz77
 * 
 */
struct lz_inst {
    byte* window = nullptr, *winStart = nullptr;
    size_t winBits = 0, winSz = 0, winPos = 0;
};

struct _match {
    hash_node<_lz_win_ref>* node;
    size_t matchSz;
};

struct match_settings {
    size_t maxLength = 0x100;
    size_t goodLength = 0xffff;
};

//gets longest match or something
_match longest_match(hash_node<_lz_win_ref>* firstMatch, lz_inst* ls, match_settings ms = {}) {
    if (!firstMatch || !ls)
        return {};

    hash_node<_lz_win_ref>* curMatch = firstMatch, *bestMatch = firstMatch;
    size_t bestSz = 0;

    byte *cur = ls->window, *comp = curMatch->val.ptr;

    //get le best match
    while (curMatch) {
        size_t mSz = 0;
        //do le matching
        do { } while (*cur++ == *comp++ && ++mSz < ms.maxLength);

        if (mSz > bestSz) {
            bestMatch = curMatch;
            bestSz = mSz;
        }

        if (mSz >= ms.goodLength) break;

        curMatch = curMatch->prev;
    }

    //return the best
    return {
        .node = bestMatch,
        .matchSz = bestSz
    };
}

//shifts lz77 window
void win_shift(lz_inst* ls) {
    if (!ls) return;
}

/**
 *
 * Main lz77 function
 *
 * lz77_encode
 *
 * idk write somethign smart for this
 * 
 * note when using function
 * 
 * THE BUFFER RETURNED, STORES DATA AS
 * 3 BYTE INTEGERS!! IN A SIZE_T* (4 byte)
 * CREATE A CLASS OR SOMETHING TO HANDLE
 * DIFFERENT SIZED STREAM LIKE THIS ONE
 * TO READ 3 BYTE VALUES INSTEAD OF 4 
 * BYTE
 * 
 * could possibly fix if you didnt do copy
 * to buffer in the micro vector but oh
 * well
 * 
 * Recommended solution: convert the size_t*
 * into a realigned buffer and then use the 
 * buffer for further parsing
 *
 */
#define _WIN_SHIFT(ls, curByte, winShift) (ls).window += (winShift); \
    (ls).winPos += (winShift); \
    (enc_byte) += (winShift)  

lzr lz77_encode(byte* dat, size_t sz, size_t winBits) {
    const size_t nMallocs = 0xf;

    //just some window constants
    const size_t lookAhead = 0x100;
    const size_t winMem = 0x7fff - lookAhead; //32kb - loohAhead = other window memory

    //size of lz77 alphabet
    const size_t lzAlphaSz = 0x120;
   
    //distance char counts for the distance tree
    const size_t dsb_len = sizeof(DistanceBase) / sizeof(i32), lnb_len = sizeof(LengthBase) / sizeof(i32);
    size_t* distCharCount = new size_t[dsb_len];
    ZeroMem(distCharCount, dsb_len);

    //encode result
    lzr res = {
        .encDat = nullptr,
        .encSz = 0,
        .winBits = winBits,
        .distTree = nullptr,
        .charCounts = nullptr
    };

    //construct important info and stuff
    lz_inst ls = {
        .window = dat,   //lz77 window ptr
        .winStart  = dat, //where the window begins
        .winBits = winBits,     //window bits (2 ^ b = s), where b is win bits and s is window size in bytes
        .winSz = (size_t) 1 << winBits,  //size of window in bytes
        .winPos = 0
    };

    //hash table
    //TODO: remove the std::string template??
    linked_map<_lz_win_ref, 15> hashTable;

    //micro vector creation / allocation
    const size_t minSz = 0xffff;
    const u32
        p_vChunkSz = sz / nMallocs,
        vChunkSz = p_vChunkSz < minSz ? minSz : p_vChunkSz;

    micro_vector<3> enc_dat = micro_vector<3>(0, vChunkSz);

    //
    byte* enc_byte = dat, * end = dat + sz;

    //allocate char count container thing
    res.charCounts = new u32[lzAlphaSz];
    ZeroMem(res.charCounts, lzAlphaSz);

    //main encode loop
    do {
        //store the last
        hash_node<_lz_win_ref>* last = hashTable.insert(reinterpret_cast<char*>(ls.window), 3, {
            .len = 3,
            .ptr = ls.window
        });

        size_t winShift = 1;
        
        if (last) {
            //then back ref
            _match m = longest_match(last, &ls);

            if (!m.node) {
                std::cout << "Bro how tf did you mess up this badly ;-;-;-;" << std::endl;
                _safe_free_a(distCharCount);
                enc_dat.free();
                return {};
            }

            if (m.matchSz <= 3) goto def_match; //make sure match is longer than 3 so it's worth it

            //encode back ref
            const size_t

                lenIdx = lz77_get_len_idx(m.matchSz),

                dist = (size_t)(ls.window - m.node->val.ptr), //get distance
                distIdx = lz77_get_dist_idx(dist); //get distance base

            //make sure things are valid
            if (distIdx >= dsb_len || distIdx < 0) {
                std::cout << "ZLib Error: invalid match distance!";
                _safe_free_a(distCharCount);
                enc_dat.free();
            }

            if (lenIdx >= lnb_len || lenIdx < 0) {
                std::cout << "ZLib Error: invalid match length!";
                _safe_free_a(distCharCount);
                enc_dat.free();
            }


            //get extras after quick error check
            const size_t distExtra = dist - DistanceBase[min(distIdx, dsb_len)],
                         lenExtra = m.matchSz - LengthBase[min(lenIdx, lnb_len)];

            distCharCount[distIdx]++; //increase char count yk

            //if we get back ref change shift amount
            winShift = m.matchSz;

            //add back ref val
            enc_dat.push(257 + lenIdx);

            /*
            
            Extra Info Format:

            1 byte: distance Index - [i]
            2 byte: dist extra     - [e]

            iee

            ----------------------

            2 byte: len extra
            
            */
            enc_dat.push(
                (distIdx << 16) & 0xff | (distExtra & 0xffff)
            );

            enc_dat.push(
                lenExtra & 0xffff
            );
        }

    def_match:

        //encode basic character
        enc_dat.push(*ls.window);
        res.charCounts[*ls.window]++;

        //window shifting
        _WIN_SHIFT(ls, curByte, winShift);

        if (ls.winPos >= winMem) {
            hashTable.clear();
            ls.winPos = 0;
        }
    } while (enc_byte < end);

    //construct distance tree and other stuff we need
    res.distTree = GenCanonicalTreeFromCounts<dsb_len>(distCharCount);
    
    if (!res.distTree)
        std::cout << "Something went wrong when generating distance tree!" << std::endl;

    //copy over data
    res.encSz = enc_dat.sz;
    enc_dat.copyToBuffer(res.encDat, res.encSz);

    //mem management and return
    _safe_free_a(distCharCount);
    enc_dat.free();
    return res;
}

/**
 *
 * DeflateBlock
 *
 * deflates a given block of data
 *
 */

enum deflate_block_type {
    dfb_none,
    dfb_static,
    dfb_dynamic
};

void WriteDeflateBlockToStream(BitStream* stream, bin* block_data, const size_t winBits, const i32 level, bool finalBlock = false) {
    deflate_block_type bType = level > 0 ? dfb_dynamic : dfb_static;

    //write some info
    stream->writeBit(finalBlock & 0b01);
    stream->writeBit((bType & 0b01));
    stream->writeBit((bType & 0b10) >> 1);

    //compresss
    switch (level) {

        //no compression so we just write stuff
    case 0: {

        break;
    }

          //TODO: huffman and other levels from 0-10
    case 1: break;

        //max compressoin
    case 2: {

        break;
    }
    default:
        return;
    }
}

/**
 *
 * Balloon::Deflate
 *
 * Core deflate function for zlib. Compresses all
 * the data and returns a result containing the
 * compressed data
 *
 */

balloon_result Balloon::Deflate(byte* data, size_t sz, u32 compressionLevel, const size_t winBits) {
    //quick error check
    if (data == nullptr || sz <= 0 || winBits > 15) return {};

    BitStream rStream = BitStream(0xff);

    //generate some of the fields
    byte cmf = (0x08) | (((winBits - 8) & 0b1111) << 4);
    byte flg = 0; //all the data is just going to be 0

    flg = 0x1f - (cmf * 0x100) % 0x1f; //compute the flag checksum

    //write the cmf and flag bytes
    rStream.writeValue(cmf);
    rStream.writeValue(flg);

    //deflate everything officially
    //TODO: multiple blocks of data



    return { 0 };
}

//line 861, lets see how off this comment gets