#pragma once

#include "DataPointer.h"
namespace MultigridOps
{

	template <int e>
	void Interpolate(
		DataPointer<(1 << e) - 1, (1 << e) - 1>& to,
		DataPointer<(1 << (e - 1)) - 1, (1 << (e - 1)) - 1>& from)
	{
		constexpr int wF = (1 << e) - 1;
		constexpr int hF = (1 << e) - 1;

		constexpr int wC = (1 << (e - 1)) - 1;
		constexpr int hC = (1 << (e - 1)) - 1;

#pragma omp parallel for
		for (int j = 0; j < hF; j++)
		{
			for (int i = 0; i < wF; i++)
			{
				//This version of interpolate is specifically for data pointers
				//that represent the interior points of a grid with haloed boundaries.
				//The halo isn't actually included in the data pointer, but the indexing
				//Works out really cleanly if we pretend that it is.
				int myI = i + 1;
				int myJ = j + 1;
				if (myI % 2 == 0 && myJ % 2 == 0)
				{
					int hisI = myI / 2;
					int hisJ = myJ / 2;
					hisI -= 1;
					hisJ -= 1;
					to.Set(i, j, from.Get(hisI, hisJ));
				}
				else if (myI % 2 == 0 && myJ % 2 == 1)
				{
					double val = 0;
					int hisI = myI / 2;
					int hisJA = myJ / 2;
					int hisJB = myJ / 2 + 1;
					hisI -= 1;
					hisJA -= 1;
					hisJB -= 1;
					if(hisJA >= 0)
						val += from.Get(hisI, hisJA) * 0.5;
					if (hisJB < hC)
						val += from.Get(hisI, hisJB) * 0.5;
					to.Set(i, j, val);
				}
				else if (myI % 2 == 1 && myJ % 2 == 0)
				{
					double val = 0;
					int hisIA = myI / 2;
					int hisIB = myI / 2 + 1;
					int hisJ = myJ / 2;
					hisIA -= 1;
					hisIB -= 1;
					hisJ -= 1;
					if (hisIA >= 0)
						val += from.Get(hisIA, hisJ) * 0.5;
					if (hisIB < wC)
						val += from.Get(hisIB, hisJ) * 0.5;
					to.Set(i, j, val);
				}
				else
				{
					double val = 0;
					int hisIA = myI / 2;
					int hisIB = myI / 2 + 1;
					int hisJA = myJ / 2;
					int hisJB = myJ / 2 + 1;
					hisIA -= 1;
					hisIB -= 1;
					hisJA -= 1;
					hisJB -= 1;
					if (hisIA >= 0 && hisJA >= 0)
						val += from.Get(hisIA, hisJA) * 0.25;
					if (hisIB < wC && hisJA >= 0)
						val += from.Get(hisIB, hisJA) * 0.25;
					if (hisIA >= 0&& hisJB < hC)
						val += from.Get(hisIA, hisJB) * 0.25;
					if (hisIB < wC && hisJB < hC)
						val += from.Get(hisIB, hisJB) * 0.25;
					to.Set(i, j, val);
				}
			}
		}
	}

	template <int e>
	void Interpolate(
		DataPointer<(1 << e) +1, (1 << e) + 1>& to,
		DataPointer<(1 << (e - 1)) + 1, (1 << (e - 1)) + 1>& from)
	{
		constexpr int wF = (1 << e) + 1;
		constexpr int hF = (1 << e) + 1;
#pragma omp parallel for
		for (int j = 0; j < hF; j++)
		{
			for (int i = 0; i < wF; i++)
			{
				double value = 0;
				int ci = i / 2;
				int cj = j / 2;
				if (i % 2 == 0 && j % 2 == 0)
				{
					value = from.Get(ci, cj);
					to.Set(i, j, value);
				}
				else if (i % 2 == 0)
				{
					value += from.Get(ci, cj) * 0.5;
					value += from.Get(ci, cj + 1) * 0.5;
					to.Set(i, j, value);
				}
				else if (j % 2 == 0)
				{
					value += from.Get(ci, cj) * 0.5;
					value += from.Get(ci + 1, cj) * 0.5;
					to.Set(i, j, value);
				}
				else
				{
					value += from.Get(ci, cj) * 0.25;
					value += from.Get(ci + 1, cj) * 0.25;
					value += from.Get(ci, cj + 1) * 0.25;
					value += from.Get(ci + 1, cj + 1) * 0.25;
					to.Set(i, j, value);
				}
			}
		}
	}

	template <int e>
	void Restrict(
		DataPointer<(1 << e) - 1, (1 << e) - 1>& from,
		DataPointer<(1 << (e - 1)) - 1, (1 << (e - 1)) - 1>& to)
	{
		constexpr int wC = (1 << (e - 1)) - 1;
		constexpr int hC = (1 << (e - 1)) - 1;
#pragma omp parallel for
		for (int j = 0; j < hC; j++)
		{
			for (int i = 0; i < wC; i++)
			{
				//This version of restrict is specifically for data pointers
				//That represent the interior points of a grid with haloed boundaries.
				//The halo isn't actually included in the data pointer, but the indexing
				//Works out really cleanly if we pretend that it is.
				int myI = i + 1;
				int myJ = j + 1;
				int hisI = 2 * myI;
				int hisJ = 2 * myJ;
				hisI -= 1;
				hisJ -= 1;

				double val = 0;
				val += from.Get(hisI - 1, hisJ - 1) / 16;
				val += from.Get(hisI, hisJ - 1) / 8;
				val += from.Get(hisI + 1, hisJ - 1) / 16;
				val += from.Get(hisI - 1, hisJ) / 8;
				val += from.Get(hisI, hisJ) / 4;
				val += from.Get(hisI + 1, hisJ) / 8;
				val += from.Get(hisI - 1, hisJ + 1) / 16;
				val += from.Get(hisI, hisJ + 1) / 8;
				val += from.Get(hisI + 1, hisJ + 1) / 16;
				to.Set(i, j, val);
			}
		}
	}

	template <int e>
	void RestrictSimple(
		DataPointer<(1 << e) + 1, (1 << e) + 1>& from,
		DataPointer<(1 << (e - 1)) + 1, (1 << (e - 1)) + 1>& to)
	{
		constexpr int wC = (1 << (e - 1)) + 1;
		constexpr int hC = (1 << (e - 1)) + 1;

#pragma omp parallel for
		for (int j = 0; j < hC; j++)
		{
			for (int i = 0; i < wC; i++)
			{
				to.Set(i, j, from.Get(2 * i, 2 * j));
			}
		}
	}

	template <int e>
	void Restrict(
		DataPointer<(1 << e) + 1, (1 << e) + 1>& from, 
		DataPointer<(1 << (e - 1)) + 1, (1 << (e - 1)) + 1>& to)
	{
		constexpr int wC = (1 << (e - 1)) + 1;
		constexpr int hC = (1 << (e - 1)) + 1;

#pragma omp parallel for
		for (int j = 0; j < hC; j++)
		{
			for (int i = 0; i < wC; i++)
			{
				if (i == 0 && j == 0)
				{
					to.Set(i, j, from.Get(2 * i, 2 * j));
				}
				else if (i == 0 && j == hC - 1)
				{
					to.Set(i, j, from.Get(2 * i, 2 * j));
				}
				else if (i == wC - 1 && j == 0)
				{
					to.Set(i, j, from.Get(2 * i, 2 * j));
				}
				else if (i == wC - 1 && j == hC - 1)
				{
					to.Set(i, j, from.Get(2 * i, 2 * j));
				}
				else if (i == 0)
				{
					double val = 0;
					val += from.Get(2 * i, 2 * j - 1) / 4;
					val += from.Get(2 * i, 2 * j) / 2;
					val += from.Get(2 * i, 2 * j + 1) / 4;
					to.Set(i, j, val);
				}
				else if (i == wC - 1)
				{
					double val = 0;
					val += from.Get(2 * i, 2 * j - 1) / 4;
					val += from.Get(2 * i, 2 * j) / 2;
					val += from.Get(2 * i, 2 * j + 1) / 4;
					to.Set(i, j, val);
				}
				else if (j == 0)
				{
					double val = 0;
					val += from.Get(2 * i - 1, 2 * j) / 4;
					val += from.Get(2 * i, 2 * j) / 2;
					val += from.Get(2 * i + 1, 2 * j) / 4;
					to.Set(i, j, val);
				}
				else if (j == hC - 1)
				{
					double val = 0;
					val += from.Get(2 * i - 1, 2 * j) / 4;
					val += from.Get(2 * i, 2 * j) / 2;
					val += from.Get(2 * i + 1, 2 * j) / 4;
					to.Set(i, j, val);
				}
				else
				{
					double val = 0;
					val += from.Get(2 * i - 1, 2 * j - 1) / 16;
					val += from.Get(2 * i, 2 * j - 1) / 8;
					val += from.Get(2 * i + 1, 2 * j - 1) / 16;
					val += from.Get(2 * i - 1, 2 * j) / 8;
					val += from.Get(2 * i, 2 * j) / 4;
					val += from.Get(2 * i + 1, 2 * j) / 8;
					val += from.Get(2 * i - 1, 2 * j + 1) / 16;
					val += from.Get(2 * i, 2 * j + 1) / 8;
					val += from.Get(2 * i + 1, 2 * j + 1) / 16;
					to.Set(i, j, val);
				}
			}
		}
	}
}