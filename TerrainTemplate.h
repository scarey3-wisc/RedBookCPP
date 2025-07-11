#pragma once

#include <cstdint>
#include <string>

class TerrainTemplate
{
	friend class TerrainType;
public:
	TerrainTemplate() : value(0), mask(0) { }
	TerrainTemplate(uint16_t s);
	TerrainTemplate(uint8_t v) : value(v), mask(0xFF) { }
	TerrainTemplate(uint8_t v, uint8_t m);
	TerrainTemplate(const TerrainTemplate& t, uint8_t v, uint8_t m);
	TerrainTemplate(const TerrainTemplate& t, const TerrainTemplate& s) : TerrainTemplate(t, s.value, s.mask) { }
	TerrainTemplate(const TerrainTemplate& t) : value(t.value), mask(t.mask) { }

	bool operator==(const TerrainTemplate& t) { return t.value == value && t.mask == mask;  }
	uint16_t BitsToShort() { return (uint16_t) ((value << 8) | (mask & 0xFF)); }
	std::string BitsToString();

	bool QualitySet(uint8_t comp) { return (mask & comp) == comp; }
	uint8_t QueryQuality(uint8_t comp) { return (uint8_t) (value & comp); }
	bool HasQuality(uint8_t comp) { return (value & comp) > 0; }
	bool MatchesQualities(uint8_t comp, uint8_t mask);
	bool IsTerrainOfType(TerrainTemplate& t) { return MatchesQualities(t.value, t.mask); }
	bool IsTerraIncognita() { return mask == 0; }

protected:
	//indicates which qualities we possess
	uint8_t value;
	//indicates which qualities have been set
	uint8_t mask;

public:

	inline static constexpr uint8_t river = 0x01; 		//0000 0001
	inline static constexpr uint8_t lake = 0x02; 		//0000 0010
	inline static constexpr uint8_t cold = 0x04; 		//0000 0100
	inline static constexpr uint8_t extreme = 0x08;		//0000 1000
	inline static constexpr uint8_t elevated = 0x10;	//0001 0000
	inline static constexpr uint8_t rough = 0x20;		//0010 0000
	inline static constexpr uint8_t dry = 0x40;			//0100 0000
	inline static constexpr uint8_t barren = 0x80;		//1000 0000

	static const TerrainTemplate MOUNTAINS;
	static const TerrainTemplate PEAKS;
	static const TerrainTemplate FLAT;
	static const TerrainTemplate ROUGH;
	static const TerrainTemplate HOT;
	static const TerrainTemplate FRIGID;
	static const TerrainTemplate WARM;
	static const TerrainTemplate COOL;
	static const TerrainTemplate DESERT;
	static const TerrainTemplate JUNGLE;
	static const TerrainTemplate PLAINS;
	static const TerrainTemplate GRASSLAND;
	static const TerrainTemplate FOREST;
	static const TerrainTemplate HILLS;
	static const TerrainTemplate TUNDRA;
	static const TerrainTemplate SNOW;
	static const TerrainTemplate LAKE;
	static const TerrainTemplate RIVER;
	static const TerrainTemplate OCEAN;
};

