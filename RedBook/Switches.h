#pragma once
namespace Switches
{
	inline constexpr double LAPLACE_EROSION_ROCK_CONSTANT = 0.001;
	inline constexpr double LAPLACE_EROSION_DEPOSITION_CONSTANT = 0.02;
	inline constexpr double LAPLACE_EROSION_SEDIMENT_CONSTANT = 0.05;

	extern bool SIMPLE_TERRAIN_TYPE_RENDERING;
	extern bool OUTLINE_MAPS;
	extern bool PAINT_VORONOI_CENTERS;

	inline constexpr bool POISSON_DENSE = false;
	inline constexpr bool POISSON_BIASED = false;
	inline constexpr bool INITIAL_EXPAND_RENDERED = false;
	inline constexpr bool INITIAL_EXPAND_EMPTY = true;
	inline constexpr bool USE_BLURRING = true;
	inline constexpr bool CLEAR_IMAGE_CACHES = false;
	inline constexpr bool PARALLEL_RENDERING = true;
	inline constexpr bool PARALLEL_CONTINENT_GEN = true;
	inline constexpr bool PHOTOGRAPH_PARALLEL_TILED_RENDER = false;

	inline constexpr int MAX_SAMPLE_POINTS_IN_LAKE = 8000;
	inline constexpr int MIN_LOCAL_MAP_PIXELS_IN_LAKE = 400;

	enum TERRAIN_GEN_ALGO
	{
		COASTAL_DISTANCE,
		COASTAL_DISTANCE_TYPE_CONSTRAINED,
		TECTONIC_UPLIFT,
		FUSED_TECTONIC_AND_DISTANCE
	};
	enum FLOW_MODEL
	{
		D4_Random,
		D8_Random,
		D_Infinity
	};
	enum PAINT_TYPE
	{
		ELEVATION_CURR,
		TERRAIN,
		TERRAIN_PRETTY,
		PHOTOGRAPHY,
		CONTOUR,
		PERLIN_DISPLAY,
		ELEV_GRADIENT,
		VORONOI_PURE,
		VORONOI_INTERPOLATED,
		VORONOI_TRIANGLES
	};
	extern PAINT_TYPE CURR_PAINT_TYPE; 
	inline constexpr FLOW_MODEL CURR_FLOW_MODEL = FLOW_MODEL::D4_Random;
	inline constexpr TERRAIN_GEN_ALGO CURR_TERRAIN_ALGO = TERRAIN_GEN_ALGO::FUSED_TECTONIC_AND_DISTANCE;

}