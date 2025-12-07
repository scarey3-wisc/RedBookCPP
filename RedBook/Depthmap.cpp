#include "Depthmap.h"
#include "WorldMap.h"
#include "RegionalMap.h"
#include "Timer.h"
#include <iostream>
#include <cmath>
#include <omp.h>

bool Depthmap::USE_OPENMP_FOR_ALLOCATION = true;


float
Depthmap::GetDepth(int px, int py)
{
	return Get(px, py);
}

void
Depthmap::SetDepth(int px, int py, float depth)
{
	Set(depth, px, py);
}

void
Depthmap::OnAllocation(DepthmapBase handle, WorldMap* world)
{
	RegionalMap* parent = world->GetRegion(
		handle.myID.regionX - RegionalMap::ORIGIN_OFFSET,
		handle.myID.regionY - RegionalMap::ORIGIN_OFFSET);
	if (parent == nullptr)
	{
		std::cout << "Allocation failed; invalid regional coordinates!" << std::endl;
		return;
	}
	if (!handle.Valid())
	{
		std::cout << "Allocation failed; invalid handle!" << std::endl;
		return;
	}
	if (handle.myID.LOD == 1)
	{
		handle.SetAll(0);
		return;
	}
	RegionalDataLoc parentLoc = handle.myID.GetLowerLOD();
	if (!world->depthmaps.DataAvailable(parentLoc, false))
	{
		handle.SetAll(0);
		return;
	}
	DepthmapBase parentHandle = world->depthmaps.GetRaster(parentLoc);
	if (!parentHandle.OkayToUse())
	{
		handle.SetAll(0);
		return;
	}

	int iDel = 0;
	int jDel = 0;
	if (handle.myID.x % 2 == 1)
		iDel = (DIM - 1) / 2;
	if (handle.myID.y % 2 == 1)
		jDel = (DIM - 1) / 2;
	//We interpolate from our parent! Be careful about the boundary / skirt.
	//I found the following reference useful:
	/*
	* N 0 1 2 3 4 5 6 7 8 N
	  0 1 2 3 4 5 6 7 8 9 A
	  S B X X X X X X X B S
	   SBXXXXXXXBS
	   0123456789A
	   N012345678N
	*/
#pragma omp parallel for if(USE_OPENMP_FOR_ALLOCATION)
	for (int index = 0; index < (DIM + 2) * (DIM + 2); index++)
	{
		int i = index % (DIM + 2);
		int j = index / (DIM + 2);
		if (i == 0 || i == DIM + 1 || j == 0 || j == DIM + 1)
		{
			//We're not worrying about the skirt right now
			handle.Set(0, index);
			continue;
		}
		int iIn = i - 1;
		int jIn = j - 1;
		if (iIn % 2 == 0 && jIn % 2 == 0)
		{
			int iPIn = iIn / 2;
			int jPIn = jIn / 2;
			int iP = iPIn + 1;
			int jP = jPIn + 1;
			double val = 0;
			val += parentHandle.Get(iP + iDel, jP + jDel);
			handle.Set((float) val, index);
		}
		else if (iIn % 2 == 1 && jIn % 2 == 0)
		{
			int iPInA = iIn / 2;
			int iPInB = iIn / 2 + 1;
			int jPIn = jIn / 2;
			int iPA = iPInA + 1;
			int iPB = iPInB + 1;
			int jP = jPIn + 1;
			double val = 0;
			val += parentHandle.Get(iPA + iDel, jP + jDel) * 0.5;
			val += parentHandle.Get(iPB + iDel, jP + jDel) * 0.5;
			handle.Set((float) val, index);
		}
		else if (iIn % 2 == 0 && jIn % 2 == 1)
		{
			int iPIn = iIn / 2;
			int jPInA = jIn / 2;
			int jPInB = jIn / 2;
			int iP = iPIn + 1;
			int jPA = jPInA + 1;
			int jPB = jPInB + 1;
			double val = 0;
			val += parentHandle.Get(iP + iDel, jPA + jDel) * 0.5;
			val += parentHandle.Get(iP + iDel, jPB + jDel) * 0.5;
			handle.Set((float) val, index);
		}
		else
		{
			int iPInA = iIn / 2;
			int iPInB = iIn / 2 + 1;
			int jPInA = jIn / 2;
			int jPInB = jIn / 2;
			int iPA = iPInA + 1;
			int iPB = iPInB + 1;
			int jPA = jPInA + 1;
			int jPB = jPInB + 1;
			double val = 0;
			val += parentHandle.Get(iPA + iDel, jPA + jDel) * 0.25;
			val += parentHandle.Get(iPA + iDel, jPB + jDel) * 0.25;
			val += parentHandle.Get(iPB + iDel, jPA + jDel) * 0.25;
			val += parentHandle.Get(iPB + iDel, jPB + jDel) * 0.25;
			handle.Set((float) val, index);
		}
	}
}
void
Depthmap::OnDeallocation(DepthmapBase handle, WorldMap* world)
{
	std::cout << "Wait, we're supposed to deallocate!" << handle.myID.regionX << std::endl;
}