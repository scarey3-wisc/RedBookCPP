#pragma once
#include "DataRasterManager.h"
#include "RegionalDataLoc.h"

inline static constexpr int HeightmapDimension = 256;

class Heightmap;

using HeightmapBase = DataRasterManager<int32_t, HeightmapDimension, RegionalDataLoc, RegionalDataLocHasher, Heightmap>::LimitedDataRasterHandle;

class Heightmap : public DataRasterManager<int32_t, HeightmapDimension, RegionalDataLoc, RegionalDataLocHasher, Heightmap>::DataRasterHandle
{
public:
	Heightmap(const DataRasterHandle& h) : DataRasterHandle(h) {}
	Heightmap(const Heightmap& h) = default;
	Heightmap& operator=(const Heightmap& h) = default;

	double GetElevation(int px, int py);
	void SetElevation(int px, int py, double height);

public:
	static bool ThreadOnAllocation() { return true; }
	static bool ThreadOnDeallocation() { return true; }

	static void OnAllocation(HeightmapBase handle, WorldMap* world);
	static void OnDeallocation(HeightmapBase handle, WorldMap* world);

	inline static constexpr int DIM = HeightmapDimension;
	static bool USE_OPENMP_FOR_ALLOCATION;
};

class HeightmapManager : public DataRasterManager<int32_t, HeightmapDimension, RegionalDataLoc, RegionalDataLocHasher, Heightmap>
{
public:
	HeightmapManager (int maxCached, WorldMap* myWorld) : DataRasterManager(maxCached, myWorld) {}
	Heightmap GetRaster(const RegionalDataLoc& id)
	{
		return Heightmap(DataRasterManager::GetRaster(id));
	}
};