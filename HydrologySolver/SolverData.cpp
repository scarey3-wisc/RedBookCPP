#include "SolverData.h"
#include <iostream>

//-----------------------------------------------------------------------------

template<int e>
void 
MultigridMesh<e>::RestrictRHSTo(MultigridMesh<e - 1>& coarser)
{
	for (int j = 0; j < coarser.h; j++)
	{
		for (int i = 0; i < coarser.w; i++)
		{
			double val = 0;
			if (i == 0 && j == 0)
			{
				val = F.Get(i * 2, j * 2);
			}
			else if (i == 0 && j == coarser.h - 1)
			{
				val = F.Get(i * 2, j * 2);
			}
			else if (i == coarser.w - 1 && j == 0)
			{
				val = F.Get(i * 2, j * 2);
			}
			else if (i == coarser.w - 1 && j == coarser.h - 1)
			{
				val = F.Get(i * 2, j * 2);
			}
			else if (i == 0)
			{
				val += F.Get(i * 2, j * 2) / 2;
				val += F.Get(i * 2, j * 2 + 1) / 4;
				val += F.Get(i * 2, j * 2 - 1) / 4;
			}
			else if (i == coarser.w - 1)
			{
				val += F.Get(i * 2, j * 2) / 2;
				val += F.Get(i * 2, j * 2 + 1) / 4;
				val += F.Get(i * 2, j * 2 - 1) / 4;
			}
			else if (j == 0)
			{
				val += F.Get(i * 2, j * 2) / 2;
				val += F.Get(i * 2 + 1, j * 2) / 4;
				val += F.Get(i * 2 - 1, j * 2) / 4;
			}
			else if (j == coarser.h - 1)
			{
				val += F.Get(i * 2, j * 2) / 2;
				val += F.Get(i * 2 + 1, j * 2) / 4;
				val += F.Get(i * 2 - 1, j * 2) / 4;
			}
			else
			{
				val += F.Get(i * 2, j * 2) / 4;
				val += F.Get(i * 2 + 1, j * 2) / 8;
				val += F.Get(i * 2 - 1, j * 2) / 8;
				val += F.Get(i * 2, j * 2 + 1) / 8;
				val += F.Get(i * 2, j * 2 - 1) / 8;
				val += F.Get(i * 2 + 1, j * 2 + 1) / 16;
				val += F.Get(i * 2 - 1, j * 2 + 1) / 16;
				val += F.Get(i * 2 + 1, j * 2 - 1) / 16;
				val += F.Get(i * 2 - 1, j * 2 - 1) / 16;
			}
			coarser.F.Set(i, j, val);
		}
	}
}

//-----------------------------------------------------------------------------

template<int e>
void 
MultigridMesh<e>::InterpolateUTo(MultigridMesh<e + 1>& finer)
{
	for (int j = 0; j < finer.h; j++)
	{
		for (int i = 0; i < finer.w; i++)
		{
			double val = 0;
			if (i % 2 == 0 && j % 2 == 0)
			{
				val = U.Get(i / 2, j / 2);
			}
			else if (i % 2 == 0)
			{
				val += U.Get(i / 2, j / 2) / 2;
				val += U.Get(i / 2, j / 2 + 1) / 2;
			}
			else if (j % 2 == 0)
			{
				val += U.Get(i / 2, j / 2) / 2;
				val += U.Get(i / 2 + 1, j / 2) / 2;
			}
			else
			{
				val += U.Get(i / 2, j / 2) / 4;
				val += U.Get(i / 2 + 1, j / 2) / 4;
				val += U.Get(i / 2, j / 2 + 1) / 4;
				val += U.Get(i / 2 + 1, j / 2 + 1) / 4;
			}
			finer.U.Set(i, j, val);
		}
	}
}