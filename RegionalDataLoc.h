#pragma once
#include <bit>
#include "HashUtils.h"
struct RegionalDataLoc
{
	uint16_t regionX, regionY; // Which RegionalMap is this heightmap part of?
	//Note that we use WORLD coordinates for this, not array coordinates.
	uint8_t x, y; // Position of this heightmap in the RegionalMap grid
	uint8_t LOD; // How many of these heightmaps fit on a single dimension of the RegionalMap? Expect 1, 2, 4, 8, or 16

	RegionalDataLoc(uint16_t rX, uint16_t rY, uint8_t x, uint8_t y, uint8_t LOD) :
		regionX(rX), regionY(rY), x(x), y(y), LOD(LOD) {
	}
	bool operator==(const RegionalDataLoc& other) const
	{
		return regionX == other.regionX &&
			regionY == other.regionY &&
			x == other.x &&
			y == other.y &&
			LOD == other.LOD;
	}
};
struct RegionalDataLocHasher
{
	uint64_t operator()(const RegionalDataLoc& k) const
	{
		uint64_t result =
			(static_cast<uint64_t>(k.regionX) << 40) |
			(static_cast<uint64_t>(k.regionY) << 24) |
			(static_cast<uint64_t>(k.x) << 16) |
			(static_cast<uint64_t>(k.y) << 8) |
			(static_cast<uint64_t>(k.LOD));
		return distribute(result);

		/*
		uint64_t seed = k.LOD;
		uint64_t x0 = k.regionX;
		uint64_t y0 = k.regionY;
		uint64_t x1 = k.x;
		uint64_t y1 = k.y;
		seed = distribute(seed);
		seed = std::rotl(seed, std::numeric_limits<uint64_t>::digits / 3) ^ distribute(x0);
		seed = std::rotl(seed, std::numeric_limits<uint64_t>::digits / 3) ^ distribute(x1);
		seed = std::rotl(seed, std::numeric_limits<uint64_t>::digits / 3) ^ distribute(y0);
		seed = std::rotl(seed, std::numeric_limits<uint64_t>::digits / 3) ^ distribute(y1);
		return seed;
		*/

	}
};