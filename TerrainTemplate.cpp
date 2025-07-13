#include "TerrainTemplate.h"

using namespace std;

TerrainTemplate::TerrainTemplate(uint16_t s)
{
	mask = (uint8_t)(s & 0xff);
	value = (uint8_t)((s >> 8) & 0xff);
}

//------------------------------------------------------------------------------

TerrainTemplate::TerrainTemplate(uint8_t v, uint8_t m)
{
	value = v; 
	mask = (uint8_t)(m | v);
}

//------------------------------------------------------------------------------

TerrainTemplate::TerrainTemplate(const TerrainTemplate& t, uint8_t v, uint8_t m)
{
	m = (uint8_t)(m | v);
	uint8_t independantMask = (uint8_t) (m ^ t.mask);
	uint8_t overlappingMask = (uint8_t) (m & t.mask);
	uint8_t valueAgreementPoints = (uint8_t) ~(v ^ t.value);
	mask = (uint8_t) (independantMask | (overlappingMask & valueAgreementPoints));
	value = (uint8_t) ((v | t.value) & mask);
}

//------------------------------------------------------------------------------

string 
TerrainTemplate::BitsToString()
{
	bool valueBits[8];
	bool maskBits[8];
	uint8_t v = value;
	uint8_t m = mask;
	for (int i = 0; i < 8; i++)
	{
		valueBits[i] = (v & 1) != 0;
		maskBits[i] = (m & 1) != 0;
		v = (uint8_t)(v >> 1);
		m = (uint8_t)(m >> 1);
	}
	string result;
	for (int i = 7; i >= 0; i--)
	{
		if (valueBits[i])
			result += "1";
		else
			result += "0";
		if (i == 4)
			result += " ";
	}
	result += " [";
	for (int i = 7; i >= 0; i--)
	{
		if (maskBits[i])
			result += "1";
		else
			result += "0";
		if (i == 4)
			result += " ";
	}
	result += "]";
	return result;
}

//------------------------------------------------------------------------------

bool
TerrainTemplate::MatchesQualities(uint8_t comp, uint8_t mask)
{
	//If I'm interpreting past Stephen correctly, the idea here is:
	//1. every bit I've marked as 'set' (ie the mask), the test must also have 'set'
	//2. for every bit I've marked as 'set', the test and I must have the same value
	if ((mask & this->mask) != mask)
		return false;
	return (comp & mask) == (value & mask);
}

//------------------------------------------------------------------------------


const TerrainTemplate TerrainTemplate::MOUNTAINS(elevated | rough, 0);
const TerrainTemplate TerrainTemplate::PEAKS(MOUNTAINS, barren, 0);
const TerrainTemplate TerrainTemplate::FLAT(0, rough);
const TerrainTemplate TerrainTemplate::ROUGH(rough, 0);
const TerrainTemplate TerrainTemplate::HOT(extreme, cold);
const TerrainTemplate TerrainTemplate::FRIGID(cold | extreme, 0);
const TerrainTemplate TerrainTemplate::WARM(0, cold | extreme);
const TerrainTemplate TerrainTemplate::COOL(cold, extreme);
const TerrainTemplate TerrainTemplate::DESERT(dry | barren | extreme, 0);
const TerrainTemplate TerrainTemplate::JUNGLE(HOT, 0, dry | barren | elevated);
const TerrainTemplate TerrainTemplate::PLAINS(dry, barren | cold);
const TerrainTemplate TerrainTemplate::GRASSLAND(WARM, 0, dry | barren);
const TerrainTemplate TerrainTemplate::FOREST(cold, dry | barren);
const TerrainTemplate TerrainTemplate::HILLS(rough, elevated);
const TerrainTemplate TerrainTemplate::TUNDRA(FRIGID, 0, barren);
const TerrainTemplate TerrainTemplate::SNOW(FRIGID, barren, 0);
const TerrainTemplate TerrainTemplate::LAKE(lake, 0);
const TerrainTemplate TerrainTemplate::RIVER(river, 0);
const TerrainTemplate TerrainTemplate::OCEAN(barren, elevated | rough | dry | extreme);
