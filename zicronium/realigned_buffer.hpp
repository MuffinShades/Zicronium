#pragma once

/**
 * 
 * Realigned Buffer
 * 
 * This class is used to take buffers that store
 * values in strange sizes (ie 3 bytes) in normal
 * sizes like 4 bytes and allow you to extract the
 * odd values easily
 * 
 * Smaller version of ByteStream that works on 
 * buffers which store values of sizes different 
 * than the C buffer type
 * 
 * IE u32 buffer holding 3 byte values
 * 
 */

#include <bit>

template<class _Ty, size_t it_sz> class realigned_buffer {
private:
	_Ty* int_buf = nullptr;
	size_t int_sz = 0;
	size_t cur = 0;
	bool lil = false;
public:
	size_t length = 0;
	realigned_buffer(_Ty* buf, size_t len) {
		this->int_buf = buf;
		this->int_sz = len;
		this->cur = 0;

		static_assert(sizeof(_Ty) > 0);

		this->computeLength();

		lil = std::endian::native == std::endian::little;
	}

	void computeLength() {
		this->length = (size_t)(this->int_sz* sizeof(_Ty)) / it_sz;
	}

	u64& operator[] (size_t idx) {
		this->computeLength();
		static_assert(idx < this->length);

		i32 tz = sizeof(_Ty), 
			sIdx = (idx * it_sz) / this->int_sz,
			subStart = (idx * it_sz) % this->int_sz,
			lSz = tz - subStart,
			rSz = this->int_sz - lSz;

		static_assert(subStart < this->int_sz - 1);

		const _Ty l = this->int_buf[subStart], r = this->int_buf[subStart+1];

		//big endian
		const u64 L = (
			l & (
				(1 << (lSz * 8)) - 1)
			),
			ts = tz * 8,
			tm = (1 << ts) - 1;

		//middle copy
		while (rSz > tz) {
			L <<= ts;
			if (subStart < this->int_sz - 1)
				L |= this->int_buf[subStart++];
			rSz -= tz;
		}

		static_assert(subStart < this->int_sz);

		//pre calculations for right side
		const _Ty Rb = this->int_buf[subStart], Rbi = rSz * 8;

		//combine left and right sides to get Big endian Value
		const u64 BeV = (L << (Rbi * 8)) |
			(
				(Rb >>
					((tz - rSz) * 8)
					)
				&
				((1 << Rbi) - 1)
				);

		if (lil)
			return BeV; // TODO: create functino to reverse bytes in a value
		else
			return BeV;
	}
};