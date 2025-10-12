#pragma once
#include "RegionalDataRasterManager.h"

class Heightmap;

using HeightmapBase = RegionalDataRasterManager<int32_t, Heightmap>::LimitedDataRasterHandle;

class Heightmap : public RegionalDataRasterManager<int32_t, Heightmap>::DataRasterHandle
{
public:
	Heightmap(const DataRasterHandle& h) : DataRasterHandle(h) {}
	Heightmap(const Heightmap& h) = default;
	Heightmap& operator=(const Heightmap& h) = default;

	double GetElevation(int px, int py);
	void SetElevation(int px, int py, double height);

public:
	static bool ThreadForAllocation() { return true; }

	static void OnAllocation(HeightmapBase handle, WorldMap* world);
	static void OnDeallocation(HeightmapBase handle, WorldMap* world);

	inline static constexpr int DIM = REGIONAL_DATA_DIMENSION;
	static bool USE_OPENMP_FOR_ALLOCATION;
};

class HeightmapManager : public MultiLODRegionalDataRasterManager<int32_t, Heightmap>
{
public:
	HeightmapManager (const std::array<int, REGIONAL_MAP_NUM_LODS>& cacheSizes, WorldMap* myWorld, RegionalDataRequestStack* pool) : 
		MultiLODRegionalDataRasterManager(cacheSizes, myWorld, pool) {}
	Heightmap GetRaster(const RegionalDataLoc& id)
	{
		return Heightmap(MultiLODRegionalDataRasterManager::GetRaster(id));
	}
};