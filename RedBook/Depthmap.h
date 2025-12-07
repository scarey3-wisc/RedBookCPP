#pragma once
#include "RegionalDataRasterManager.h"

class Depthmap;

using DepthmapBase = RegionalSkirtedCellDataRasterManager<float, Depthmap>::LimitedDataRasterHandle;

class Depthmap : public RegionalSkirtedCellDataRasterManager<float, Depthmap>::DataRasterHandle
{
public:
	Depthmap(const DataRasterHandle& h) : DataRasterHandle(h) {}
	Depthmap(const Depthmap& h) = default;
	Depthmap& operator=(const Depthmap& h) = default;

	float GetDepth(int px, int py);
	void SetDepth(int px, int py, float height);

public:
	static bool ThreadForAllocation() { return true; }

	static void OnAllocation(DepthmapBase handle, WorldMap* world);
	static void OnDeallocation(DepthmapBase handle, WorldMap* world);

	inline static constexpr int DIM = REGIONAL_DATA_DIM;
	static bool USE_OPENMP_FOR_ALLOCATION;
};

class DepthmapManager : public MultiLODRegionalSkirtedCellDataRasterManager<float, Depthmap>
{
public:
	DepthmapManager(const std::array<int, REGIONAL_MAP_NUM_LODS>& cacheSizes, WorldMap* myWorld, RegionalDataRequestStack* pool) :
		MultiLODRegionalDataRasterManager(cacheSizes, myWorld, pool) {
	}

	//The difference between get and demand is that Get will return an unusable heightmap
	//if a usable one isn't ready; demand will block until a usable one is ready.
	Depthmap GetRaster(const RegionalDataLoc& id)
	{
		return Depthmap(MultiLODRegionalDataRasterManager::GetRaster(id));
	}
	Depthmap DemandRaster(const RegionalDataLoc& id)
	{
		return Depthmap(MultiLODRegionalDataRasterManager::DemandRaster(id));
	}
};