#include <iostream>
#include "balloon.hpp"
#include "micro_vector.h"
#include "realigned_buffer.hpp"

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

	return 0;
}