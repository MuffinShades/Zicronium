#pragma once
#include <iostream>
#include <cstring>
#include <vector>

typedef char i8;
typedef unsigned char u8;
typedef unsigned char byte;

typedef int16_t i16;
typedef uint16_t u16;

typedef int32_t i32;
typedef uint32_t u32;

typedef int64_t i64;
typedef uint64_t u64;

template<class _Ty> inline static void ZeroMem(_Ty* dat, size_t sz) {
    if (dat == nullptr) return;
    memset(dat, 0, sizeof(_Ty) * sz);
}

inline static void ZeroMem(void* dat, size_t sz) {
    if (dat == nullptr) return;
    memset(dat, 0, sz);
}

inline static size_t GetNumSz(const unsigned long long num) {
    unsigned long long compare = 1;
    size_t nBytes = 0;
    const size_t mxSz = sizeof(unsigned long long);

    while (num >= compare) {
        nBytes++;
        compare <<= 8;
    }

    return nBytes;
}

inline static float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

struct Point {
    float x, y;
};

#define GMask(l) ((1 << (l)) - 1)

static inline std::vector<std::string> SplitString(std::string str, const char delim) {
    std::vector<std::string> res;
    const size_t len = str.length();
    std::string collector = "";
    const char* cStr = str.c_str();

    for (size_t i = 0; i < len; i++) {
        if (cStr[i] == delim) {
            res.push_back(collector);
            collector = "";
        }
        else
            collector += cStr[i];
    }

    res.push_back(collector);

    return res;
}

#define MAKE_LONG_BE(b0, b1, b2, b3, b4, b5, b6, b7) \
    ((b0) << 56) | \
    ((b1) << 48) | \
    ((b2) << 40) | \
    ((b3) << 32) | \
    ((b4) << 24) | \
    ((b5) << 16) | \
    ((b6) <<  8) | \
    ((b7) <<  0)

#define MAKE_LONG_LE(b0, b1, b2, b3, b4, b5, b6, b7) \
    ((b7) << 56) | \
    ((b6) << 48) | \
    ((b5) << 40) | \
    ((b4) << 32) | \
    ((b3) << 24) | \
    ((b2) << 16) | \
    ((b1) <<  8) | \
    ((b0) <<  0)

#define MAKE_INT_BE(b0, b1, b2, b3) \
    ((b0) << 24) | \
    ((b1) << 16) | \
    ((b2) <<  8) | \
    ((b3) <<  0)

#define MAKE_INT_LE(b0, b1, b2, b3) \
    ((b3) << 24) | \
    ((b2) << 16) | \
    ((b1) <<  8) | \
    ((b0) <<  0)

#define MAKE_SHORT_BE(b0, b1) \
    ((b0) <<  8) | \
    ((b1) <<  0)

#define MAKE_SHORT_LE(b0, b1) \
    ((b1) <<  8) | \
    ((b0) <<  0)

static inline u64 modifyByte(u64 val, i32 b, byte v) {
    u64 mask = ((1 << (sizeof(u64) << 3)) - 1);
    b <<= 3;
    mask ^= 0xff << b;
    val &= mask;
    return val | (v << b);
}

static bool _strCompare(std::string s1, std::string s2, bool lcmp = true, size_t slen = 0) {
    //length comparision
    if (lcmp)
        if (s1.length() != s2.length())
            return false;

    //compare chars
    const char* sc1 = s1.c_str(), * sc2 = s2.c_str();
    size_t sl = slen > 0 ? slen : s1.length();

    while (sl--)
        if (*sc1++ != *sc2++)
            return false;

    return true;
}

template<typename _Ty> static inline void _safe_free_a(_Ty* m) {
    if (m == nullptr) return;
    else {
        try {
            delete[] m;
            m = nullptr;
        }
        catch (std::exception e) {
            //oof
            m = nullptr;
            return;
        }
    }
}

template<typename _Ty> static inline void _safe_free_b(_Ty* m) {
    if (m == nullptr) return;
    else {
        try {
            delete m;
            m = nullptr;
        }
        catch (std::exception e) {
            //oof
            m = nullptr;
            return;
        }
    }
}

template<typename _Ty> static _Ty ArrMax(_Ty* arr, size_t len) {
    auto max = 0;

    for (size_t i = 0; i < len; i++)
        if (arr[i] > max) max = arr[i];

    return max;
}

static u64 NumReverse(u64 v, size_t bSz) {
    u64 r = 0;

    for (i32 i = 0; i < bSz; i++, r <<= 8, v >>= 8)
        r |= (v & 0xff);

    return r >> 8;
}

#ifndef max
#define max(a, b) ((a) > (b)) ? (a) : (b)
#endif

#ifndef min
#define min(a, b) ((a) < (b)) ? (a) : (b)
#endif

/**
 *
 * Void buffer macros to make void buffers 
 * less annoying
 * 
 * VOID_BUF_ADD -> adds a value to the buffer
 * VOID_BUF_GET -> gets current value in a buffer
 * VOID_BUF_INC -> goes to next element in buffer
 * VOID_BUF_SET -> set current value in a void buffer
 * VOID_BUF_SI -> basically a void buf set followed by a void buf inc
 * 
 */
#define VOID_BUF_ADD(buf, val) (buf) = (char*)(buf) + 1; \
                                *((char*)(buf)) = (val)
#define VOID_BUF_GET(buf) (*((char*)(buf)))
#define VOID_BUF_INC(buf) (buf) = (char*)(buf) + 1
#define VOID_BUF_SET(buf, val) *((char*)(buf)) = (val)
#define VOID_BUF_SI(buf, val) *((char*)(buf)) = (val); \
                               (buf) = (char*)(buf) + 1