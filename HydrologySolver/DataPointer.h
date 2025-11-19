#pragma once
#include <algorithm>

template<int w, int h>
class DataPointer
{
public:

//-----------------------------------------------------------------------------

	static DataPointer AllocateNewPointer()
	{
		return DataPointer();
	}

//-----------------------------------------------------------------------------

	void Deallocate()
	{
		delete ptr;
	}

//-----------------------------------------------------------------------------

	inline double Get(int x, int y)
	{
		return ptr[y * w + x];
	}

//-----------------------------------------------------------------------------

	inline void Set(int x, int y, double v)
	{
		ptr[y * w + x] = v;
	}

//-----------------------------------------------------------------------------

	double Norm()
	{
		double result = 0;
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				result += Get(i, j) * Get(i, j);
			}
		}
		return sqrt(result);
	}

//-----------------------------------------------------------------------------

	void Set(double value)
	{
		std::fill(ptr, ptr + w * h, value);
	}

//-----------------------------------------------------------------------------

	void Scale(double scale)
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				Set(i, j, Get(i, j) * scale);
			}
		}
	}

//-----------------------------------------------------------------------------

	void Add(DataPointer& b, double scale)
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				Set(i, j, Get(i, j) + scale * b.Get(i, j));
			}
		}
	}

//-----------------------------------------------------------------------------

	void AddCorrectionToHaloedData(DataPointer<w - 2, h - 2>& b, double scale, double min, double max)
	{
		for (int j = 1; j < h - 1; j++)
		{
			for (int i = 1; i < w - 1; i++)
			{
				double result = Get(i, j) + scale * b.Get(i - 1, j - 1);
				result = std::clamp(result, min, max);
				Set(i, j, result);
			}
		}
	}
private:
	DataPointer()
	{
		ptr = new double[w * h];
	}
private:
	double* ptr;
};

