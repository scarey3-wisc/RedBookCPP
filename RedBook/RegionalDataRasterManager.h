#pragma once
#include <bit>
#include "HashUtils.h"
#include "DataRasterManager.h"
#include <array>
#include "TileSizes.h"

inline static constexpr int REGIONAL_DATA_DIM = TileSizes::PIXELS_PER_REGIONAL_DATA_RASTER;
inline static constexpr int REGIONAL_MAP_NUM_LODS = 5; //1, 2, 4, 8, 16; remember, a LOD is "how many 256 pixel heightmaps fill the regional map?



struct RegionalDataLoc
{
	uint16_t regionX, regionY; // Which RegionalMap is this heightmap part of?
	//Note that we use WORLD coordinates for this, not array coordinates.
	uint8_t x, y; // Position of this heightmap in the RegionalMap grid
	uint8_t LOD; // How many of these heightmaps fit on a single dimension of the RegionalMap? Expect 1, 2, 4, 8, or 16

	RegionalDataLoc() :
		regionX(0), regionY(0), x(0), y(0), LOD(-1) 
	{
	}
	RegionalDataLoc(uint16_t rX, uint16_t rY, uint8_t x, uint8_t y, uint8_t LOD) :
		regionX(rX), regionY(rY), x(x), y(y), LOD(LOD) 
	{
	}
	const char* ToString() const
	{
		static thread_local char buffer[100];
		snprintf(buffer, 100, "Region (%d, %d), Pos (%d, %d), LOD %d", regionX, regionY, x, y, LOD);
		return buffer;
	}
	bool operator==(const RegionalDataLoc& other) const
	{
		return regionX == other.regionX &&
			regionY == other.regionY &&
			x == other.x &&
			y == other.y &&
			LOD == other.LOD;
	}
	int GetPriority()
	{
		return LOD;
	}
	int GetNumSectionsPerSide() const
	{
		return 1 << (LOD - 1);
	}
	inline void GetWorldCoordinates(int i, int j, double& wX, double& wY)
	{
		double wcWidth = 1.0 * TileSizes::LOCAL_MAPS_PER_REGIONAL_MAP / GetNumSectionsPerSide();
		double baseWorldX = regionX * TileSizes::LOCAL_MAPS_PER_REGIONAL_MAP + wcWidth * x;
		double baseWorldY = regionY * TileSizes::LOCAL_MAPS_PER_REGIONAL_MAP + wcWidth * y;
		double pixelWidth = wcWidth / (REGIONAL_DATA_DIM - 1);
		wX = baseWorldX + (i - 1) * pixelWidth;
		wY = baseWorldY + (j - 1) * pixelWidth;
	}
	RegionalDataLoc GetLowerLOD()
	{
		if (LOD == 1)
			return RegionalDataLoc(regionX, regionY, x, y, LOD);
		return RegionalDataLoc(regionX, regionY, x / 2, y / 2, LOD - 1);
	}
	std::array<RegionalDataLoc, 4> GetHigherLODs()
	{
		return {
			RegionalDataLoc(regionX, regionY, x * 2, y * 2, LOD + 1),
			RegionalDataLoc(regionX, regionY, x * 2 + 1, y * 2, LOD + 1),
			RegionalDataLoc(regionX, regionY, x * 2, y * 2 + 1, LOD + 1),
			RegionalDataLoc(regionX, regionY, x * 2 + 1, y * 2 + 1, LOD + 1)
		};
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

using RegionalDataRequestStack = RequestStack <RegionalDataLoc, RegionalDataLocHasher>;

class RegionalDataCheck {
public:
	virtual bool DataAvailable(const RegionalDataLoc& id, bool launchRequest) = 0;

	virtual ~RegionalDataCheck() = default; // always for polymorphic base classes
};

template <typename DataType, int w, int h, typename Policy = DataRasterDefaultPolicy>
using RegionalDataRasterManager = DataRasterManager<DataType, w, h, RegionalDataLoc, RegionalDataLocHasher, Policy>;

template <typename DataType, int w, int h, typename Policy = DataRasterDefaultPolicy>
using RegionalDataHandle = RegionalDataRasterManager<DataType, w, h, Policy>::DataRasterHandle;

template <typename DataType, int w, int h, typename Policy = DataRasterDefaultPolicy>
class MultiLODRegionalDataRasterManager : public RegionalDataCheck
{

private:
	std::vector<std::unique_ptr<RegionalDataRasterManager<DataType, w, h, Policy>>> managers;

public:
	MultiLODRegionalDataRasterManager(const std::array<int, REGIONAL_MAP_NUM_LODS>& cacheSizes, WorldMap* myWorld, RegionalDataRequestStack* pool)
	{
		for (int i = 0; i < REGIONAL_MAP_NUM_LODS; i++)
		{
			managers.emplace_back(std::make_unique<RegionalDataRasterManager<DataType, w, h, Policy>>(cacheSizes[i], myWorld, pool));
		}
	}

	// Disable copying:
	MultiLODRegionalDataRasterManager(const MultiLODRegionalDataRasterManager&) = delete;
	MultiLODRegionalDataRasterManager& operator=(const MultiLODRegionalDataRasterManager&) = delete;

	RegionalDataHandle<DataType, w, h, Policy> GetRaster(const RegionalDataLoc& id)
	{
		RegionalDataLoc lower = id;
		for (int i = id.LOD - 1; i >= 1; i--)
		{
			lower = lower.GetLowerLOD();
			managers[i - 1]->RefreshRaster(lower);
		}
		return managers[id.LOD - 1]->GetRaster(id);
	}
	RegionalDataHandle<DataType, w, h, Policy> DemandRaster(const RegionalDataLoc& id)
	{
		RegionalDataLoc lower = id;
		for (int i = id.LOD - 1; i >= 1; i--)
		{
			lower = lower.GetLowerLOD();
			managers[i - 1]->RefreshRaster(lower);
		}
		return managers[id.LOD - 1]->DemandRaster(id);
	}
	bool DataAvailable(const RegionalDataLoc& id, bool launchRequest)
	{
		if (launchRequest)
		{
			RegionalDataHandle<DataType, w, h, Policy> raster = GetRaster(id);
			return raster.OkayToUse();
		}
		else
		{
			return managers[id.LOD - 1]->RasterAvailable(id);
		}
	}

	int GetMaxCapacity(int LOD) const { return managers[LOD - 1]->GetMaxCapacity(); }
	int GetNumUsed(int LOD) const { return managers[LOD - 1]->GetNumUsed(); }
};

template <typename DataType, typename Policy = DataRasterDefaultPolicy>
using RegionalSkirtedCellDataRasterManager = RegionalDataRasterManager<DataType, REGIONAL_DATA_DIM + 2, REGIONAL_DATA_DIM + 2, Policy>;

template <typename DataType, typename Policy = DataRasterDefaultPolicy>
using RegionalSkirtedCellDataHandle = RegionalSkirtedCellDataRasterManager<DataType, Policy>::DataRasterHandle;

template <typename DataType, typename Policy = DataRasterDefaultPolicy>
using MultiLODRegionalSkirtedCellDataRasterManager = MultiLODRegionalDataRasterManager<DataType, REGIONAL_DATA_DIM + 2, REGIONAL_DATA_DIM + 2, Policy>;