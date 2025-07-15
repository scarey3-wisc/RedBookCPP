#pragma once
#include <cstdint> // for uint8_t

enum class WatermapValue : uint8_t {
	NotWater = 0,
	Ocean = 1,
	Lake = 2,
	Unknown = 3
};