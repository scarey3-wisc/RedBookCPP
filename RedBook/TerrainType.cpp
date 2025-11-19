#include "TerrainType.h"

void 
TerrainType::ActivateQuality(uint8_t qua)
{
	value = (uint8_t) (value | qua);
	mask = (uint8_t) (mask | qua);
}

//------------------------------------------------------------------------------

void 
TerrainType::DeactivateQuality(uint8_t qua)
{
	value = (uint8_t) (value & ~qua);
	mask = (uint8_t) (mask | qua);
}

//------------------------------------------------------------------------------

void 
TerrainType::SetQualityUnknown(uint8_t qua)
{
	value = (uint8_t) (value & ~qua);
	mask = (uint8_t) (mask & ~qua);
}

//------------------------------------------------------------------------------

void 
TerrainType::ApplyTerrain(const TerrainTemplate& t)
{
	uint8_t setToTrue = (uint8_t)(t.mask & t.value);
	uint8_t setToFalse = (uint8_t)(t.mask & ~t.value);
	value = (uint8_t)(value | setToTrue);
	value = (uint8_t)(value & ~setToFalse);
	mask = (uint8_t)(mask | t.mask);
}