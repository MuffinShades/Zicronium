#include "linked_map.hpp"
/*
//alloc
template<class _st, size_t hBits> void linked_map<_st, hBits>::allocHash() {
    this->free();
    this->hashSz = 1 << this->hashBits;
    this->roots = new hash_node * [this->hashSz];
}

//compute le hash
template<class _st, size_t hBits> u64 linked_map<_st, hBits>::computeHash(char* dat, const size_t len) {
    const u64 mask = (1 << this->hashBits) - 1;

    u64 hash = 0;

    for (size_t i = 0; i < len; i++) {
        hash += dat[i] ^ (i ^ dat[0] + 0xf6) * 7;
        hash ^= ~mask;
    }

    return (hash % mask) & mask;
}

template<class _storeType, size_t _> hash_node<_storeType>* linked_map<_storeType, _>::_insert(hash_node<_storeType>* n) {
    n->prev = this->roots[n->hash];
    this->roots[n->hash] = n;
    return n->prev;
}

template<class _st, size_t _> void linked_map<_st, _>::free() {
    this->clear();
    if (this->roots != nullptr) {
        delete[] this->roots;
        this->roots = nullptr;
    }
}

template<class _storeType, size_t _> hash_node<_storeType>* linked_map<_storeType, _>::insert(char* key, const size_t key_sz, const _storeType dat) {
    u64 hsh = this->computeHash(key, key_sz);
    return this->_insert(new hash_node<_storeType>{
        .hash = hsh,
        .key = key,
        .sz = key_sz,
        .val = dat
        });
}

template<class _storeType, size_t _> hash_node<_storeType>* linked_map<_storeType, _>::insert(std::string key, const _storeType dat) {
    return this->insert(
        const_cast<char*>(key.c_str()),
        key.length(),
        dat
    );
}

template<class _storeType, size_t hBits> template<class _KeyTy> hash_node<_storeType>* linked_map<_storeType, hBits>::insert(_KeyTy key, const _storeType dat) {
    const size_t k_sz = sizeof(_KeyTy);
    void* k_mem = (void*)&key;

    return this->insert(
        reinterpret_cast<char*>(k_mem),
        k_sz,
        dat
    );
}

template<class _storeType, size_t hBits> void linked_map<_storeType, hBits>::clear() {
    for (size_t i = 0; i < this->hashSz; i++) {
        hash_node<_storeType>* cur = this->roots[i];
        while (cur) {
            hash_node<_storeType>* prev = cur->prev;
            delete cur;
            cur = prev;
        }

        this->roots[i] = nullptr;
    }
}

template<class _storeType, size_t hBits> linked_map<_storeType, hBits>::linked_map() {
    this->allocHash();
}*/