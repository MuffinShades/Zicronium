#include <iostream>
#include "balloon.hpp"
#include "micro_vector.h"
#include "realigned_buffer.hpp"
#include "filewrite.hpp"

i32 main() {
	std::cout << "Micro Vector Testing" << std::endl;

	micro_vector<3> test_v;

	test_v.push(5);
	test_v.push(6);
	test_v.push(7);

	u32* testBuf = new u32[test_v.sz];

	ZeroMem(testBuf, test_v.sz);

	test_v.copyToBuffer(testBuf, test_v.sz);

	std::cout << "Res: " << std::endl;

	for (size_t i = 0; i < test_v.sz; i++)
		std::cout << testBuf[i] << std::endl;


	std::cout << "Realign Test:" << std::endl;

	realigned_buffer realigendBuf = realigned_buffer<u32, 3>(testBuf, test_v.sz);


	std::cout << "Vals: " << std::endl;

	for (size_t i = 0; i < test_v.sz; i++)
		std::cout << realigendBuf[i] << std::endl;


	std::cout << "------ Deflate Testing --------" << std::endl;

	byte testData[] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0xb, 0xb, 0xc, 0xc, 0xc, 0xc, 0xc, 0xc, 0xc, 0x69};

	//balloon_result testRes = Balloon::Deflate(testData, 20);

	//std::cout << "Result Size: " << testRes.sz << std::endl;

	//foreach_ptr(byte, val, testRes.data, testRes.sz)
		//std::cout << (int)*val << " ";

	//std::cout << std::endl;

	//balloon_result inflateTest = Balloon::Inflate(testRes.data, testRes.sz);

	//std::cout << "Inflate Size: " << inflateTest.sz << std::endl;
	
	file testFile = FileWrite::readFromBin("C:\\TestStuff\\zlib2.txt");
	balloon_result testFileDeflate = Balloon::Deflate(testFile.dat, testFile.len);

	std::cout << "testFile Deflate Size: " << testFileDeflate.sz << std::endl;
	std::cout << "checksum: " << testFileDeflate.checksum << std::endl;

	FileWrite::writeToBin("C:\\TestStuff\\zlib2_deflate.bin", testFileDeflate.data, testFileDeflate.sz);

	balloon_result testFileInflate = Balloon::Inflate(testFileDeflate.data, testFileDeflate.sz);

	std::cout << " testFile Inflate Size: " << testFileInflate.sz << std::endl;
	std::cout << "inflate checksum: " << testFileInflate.checksum << std::endl;

	FileWrite::writeToBin("C:\\TestStuff\\zlib2_inflate.txt", testFileInflate.data, testFileInflate.sz);

	
	//TODO: make sure when something goes wrong everything errors out
	//ie -> rle encode error -> all functions prior error out / return
	//TODO: find a way to cap the bit length of a character when generating a tree
	//possible solution: when generating the trees if a bitlengths greater than the max is noticed then after
	//tree generation call a fixer function that basically takes all bitlengths, orders by lengths, and adds to
	//the longer codes until all bit lengths are uniform and the sum of all bitlengths is the same
	file testFile2 = FileWrite::readFromBin("C:\\TestStuff\\wood.wld");
	balloon_result testFileDeflate2 = Balloon::Deflate(testFile2.dat, testFile2.len);

	std::cout << "testFile Deflate Size: " << testFileDeflate2.sz << std::endl;
	std::cout << "checksum: " << testFileDeflate2.checksum << std::endl;

	FileWrite::writeToBin("C:\\TestStuff\\wood.bal", testFileDeflate2.data, testFileDeflate2.sz);

	balloon_result testFileInflate2 = Balloon::Inflate(testFileDeflate2.data, testFileDeflate2.sz);

	std::cout << " testFile Inflate Size: " << testFileInflate2.sz << std::endl;
	std::cout << "inflate checksum: " << testFileInflate2.checksum << std::endl;

	FileWrite::writeToBin("C:\\TestStuff\\wood2.wld", testFileInflate2.data, testFileInflate2.sz);
	
	

	return 0;
}