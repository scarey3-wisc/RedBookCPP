#pragma once  
#include <limits>  
namespace RedBook  
{  
	inline static constexpr int TILE_RENDERING_CACHES = 6;  
	inline static constexpr int TILE_RENDERING_MINS[] = {  
		0, 401, 1601, 6401, 25601, 102401  
	};  
	inline static constexpr int TILE_RENDERING_MAXES[] = {  
		400, 1600, 6400, 25600, 102400, (std::numeric_limits<int>::max)()  
	};  
	inline static constexpr long TILE_RENDERING_CAP[] = {  
		20000000l, 5000000l, 5000000l, 5000000l, 10000000l, 40000000l  
	};  
	inline static constexpr const char* TILE_RENDERING_NAMES[] = {  
		"Tiny", "Small", "Medium", "Large", "Massive", "Huge"  
	};  
}