#pragma once
#include <iostream>
#include <cstring>
#include "msutil.hpp"

#define MSFL_EXP

//bit stream class for reading and doing stuff with bits
class BitStream {
public:
	u32* bytes = nullptr;
	i32 pos = 0, rPos = 0;
	i32 lBit = 0;
	u32 cByte = 0;
	u32 rBit = 0;
	size_t sz, bsz, asz;
	//bit stream for reading existing bytes / bits
	MSFL_EXP BitStream(u32* bytes, size_t len);
	//basic zero allocation bit stream for writing
	MSFL_EXP BitStream(size_t len);
	MSFL_EXP i32 readBit();
	MSFL_EXP i32 readByte();
	MSFL_EXP i32 readNBits(i32 nBits);
	MSFL_EXP i32 readValue(size_t size);
	MSFL_EXP static u32* bitsToBytes(u32* bits, size_t len);
	MSFL_EXP static u32* bitsToBytes(u32 bits, size_t len);
	MSFL_EXP void checkWritePosition();
	MSFL_EXP void writeBit(u32 bit);
	//write multiple bits
	MSFL_EXP void writeNBits(u32* bits, size_t len);
	//write short, long, int, uint, etc.
	MSFL_EXP void writeValue(byte val);
	MSFL_EXP void writeValue(short val);
	MSFL_EXP void writeValue(i32 val);
	MSFL_EXP void writeValue(u32 val);
	MSFL_EXP void writeValue(i64 val);
	MSFL_EXP void writeValue(u64 val);
	//allocate new bit chunk
	MSFL_EXP void allocNewChunk();
	MSFL_EXP void seek(i32 pos);
	MSFL_EXP void seekw(i32 wPos);
	MSFL_EXP i32 tell();
	MSFL_EXP i32 tellw();
	//remove unused bytes
	MSFL_EXP void clip();
	//allocation function
	MSFL_EXP void calloc(size_t sz);
};

struct balloon_result {
	byte* data;
	u32 checksum;
	size_t sz;
	byte compressionMethod;
};

class Balloon {
	balloon_result Deflate(byte* data, size_t sz, u32 compressionLevel = 2, const size_t winBits = 0xf);
	balloon_result Inflate(byte* data, size_t sz);
};