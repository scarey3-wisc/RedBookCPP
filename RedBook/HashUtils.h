#pragma once
template<typename T>
T xorshift(const T& n, int i) {
	return n ^ (n >> i);
}
inline uint32_t distribute(const uint32_t& n) {
	uint32_t p = 0x55555555ul; // pattern of alternating 0 and 1
	uint32_t c = 3423571495ul; // random uneven integer constant; 
	return c * xorshift(p * xorshift(n, 16), 16);
}
inline uint64_t distribute(const uint64_t& n) {
	uint64_t p = 0x5555555555555555ull; // pattern of alternating 0 and 1
	uint64_t c = 17316035218449499591ull;// random uneven integer constant; 
	return c * xorshift(p * xorshift(n, 32), 32);
}