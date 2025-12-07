#pragma once
namespace TileSizes
{
	inline static constexpr int METERS_PER_LOCAL_MAP = 10240;
	inline static constexpr int LOCAL_MAPS_PER_REGIONAL_MAP = 1 << (5 - 1); //2^4 = 16
	inline static constexpr int METERS_PER_REGIONAL_MAP = METERS_PER_LOCAL_MAP * LOCAL_MAPS_PER_REGIONAL_MAP;
	inline static constexpr int VORONOI_GRID_CELLS_PER_REGIONAL_MAP = 192;
	inline static constexpr int REGIONAL_MAPS_OFFSET_FROM_ORIGIN = 500;
	inline static constexpr int PIXELS_PER_REGIONAL_DATA_RASTER = 257;
}