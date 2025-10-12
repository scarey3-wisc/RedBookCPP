#include "Heightmap.h"
#include "WorldMap.h"
#include "RegionalMap.h"
#include "Timer.h"
#include <iostream>
#include <cmath>
#include <omp.h>

bool Heightmap::USE_OPENMP_FOR_ALLOCATION = true;

namespace
{
	const double Shift = std::pow(2, 16);
}

double 
Heightmap::GetElevation(int px, int py)
{
	int32_t stored = Get(px, py);
	double calculated = (double)stored;
	return calculated / Shift;
}
void 
Heightmap::SetElevation(int px, int py, double height)
{
	height *= Shift;
	int32_t stored = (int)height;
	Set(stored, px, py);
}

void 
Heightmap::OnAllocation(HeightmapBase handle, WorldMap* world)
{
	RegionalMap* parent = world->GetRegion(
		handle.myID.regionX - RegionalMap::ORIGIN_OFFSET, 
		handle.myID.regionY - RegionalMap::ORIGIN_OFFSET);
	if(parent == nullptr)
	{
		std::cout << "Allocation failed; invalid regional coordinates!" << std::endl;
		return;
	}
	if (!handle.Valid())
	{
		std::cout << "Allocation failed; invalid handle!" << std::endl;
		return;
	}
	//width in world coordinates - probably between 1 and RegionalMap::DIMENSION, in powers of 2
	double wcWidth = RegionalMap::DIMENSION / handle.myID.GetNumSectionsPerSide();
	double baseWorldX = handle.myID.regionX * RegionalMap::DIMENSION + wcWidth * handle.myID.x;
	double baseWorldY = handle.myID.regionY * RegionalMap::DIMENSION + wcWidth * handle.myID.y;
	double pixelWidth = wcWidth / DIM;

	if(USE_OPENMP_FOR_ALLOCATION)
	{
#pragma omp parallel for
		for (int index = 0; index < DIM * DIM; index++)
		{
			/*if (index == 0)
			{
				std::cout << "Allocating heightmap for region " << handle.myID.regionX << ", " << handle.myID.regionY
					<< " at LOD " << (int)handle.myID.LOD << " (" << baseWorldX << ", " << baseWorldY << " to "
					<< (baseWorldX + wcWidth) << ", " << (baseWorldY + wcWidth) << ") " << " using "
					<< omp_get_num_threads() << " threads." << std::endl;
			}*/
			int i = index % DIM;
			int j = index / DIM;
			double wX = baseWorldX + i * pixelWidth;
			double wY = baseWorldY + j * pixelWidth;
			double elev = parent->CalculateElevation(wX, wY);
			elev *= Shift;
			int32_t stored = (int)elev;
			handle.Set(stored, index);
		}
	}
	else
	{
		for (int index = 0; index < DIM * DIM; index++)
		{
			int i = index % DIM;
			int j = index / DIM;
			double wX = baseWorldX + i * pixelWidth;
			double wY = baseWorldY + j * pixelWidth;
			double elev = parent->CalculateElevation(wX, wY);
			elev *= Shift;
			int32_t stored = (int)elev;
			handle.Set(stored, index);
		}
	}

}
void 
Heightmap::OnDeallocation(HeightmapBase handle, WorldMap* world)
{
	std::cout << "Wait, we're supposed to deallocate!" << handle.myID.regionX << std::endl;
}