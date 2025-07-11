#pragma once
#include "TerrainTemplate.h"
class TerrainType : public TerrainTemplate
{
public:
	TerrainType() : TerrainTemplate() { }
	TerrainType(uint16_t s) : TerrainTemplate(s) { }
	TerrainType(uint8_t v) : TerrainTemplate(v) { }

	TerrainType(uint8_t v, uint8_t m) : TerrainTemplate(v, m) { }

	TerrainType(const TerrainTemplate& t, uint8_t v, uint8_t m) : TerrainTemplate(t, v, m) { }
	TerrainType(const TerrainTemplate& t, const TerrainTemplate& s) : TerrainTemplate(t, s) { }
	TerrainType(const TerrainTemplate& t) : TerrainTemplate(t) { }

	void ActivateQuality(uint8_t qua);
	void DeactivateQuality(uint8_t qua);
	void SetQualityUnknown(uint8_t qua);
	void ApplyTerrain(const TerrainTemplate& t);
};

