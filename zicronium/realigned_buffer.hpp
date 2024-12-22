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

template<class _Ty, size_t it_sz> class realigned_buffer {
private:
	_Ty* int_buf = nullptr;
	size_t int_sz = 0;
	size_t cur = 0;
public:

	realigned_buffer(_Ty* buf, size_t len) {
		this->int_buf = buf;
		this->int_sz = len;
		this->cur = 0;

		static_assert(sizeof(_Ty) > 0);
	}

	u64 readVal() {

	}
};