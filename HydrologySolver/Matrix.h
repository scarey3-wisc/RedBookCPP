#pragma once

#include <algorithm>

template <int width, int height>
struct Matrix
{

public:

	//-----------------------------------------------------------------------------

	static Matrix AllocateNewPointer()
	{
		return Matrix();
	}

	//-----------------------------------------------------------------------------

	void Deallocate()
	{
		delete ptr;
	}

	//-----------------------------------------------------------------------------

	inline double Get(int x, int y) const
	{
		return ptr[y * width + x];
	}

	//-----------------------------------------------------------------------------

	inline void Set(int x, int y, double v)
	{
		ptr[y * width + x] = v;
	}

	//-----------------------------------------------------------------------------

	inline void Add(int x, int y, double v)
	{
		ptr[y * width + x] += v;
	}

	//-----------------------------------------------------------------------------

	void SetAll(double value)
	{
		std::fill(ptr, ptr + width * height, value);
	}

	//-----------------------------------------------------------------------------

	void Scale(double scale)
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				Set(i, j, Get(i, j) * scale);
			}
		}
	}

	//-----------------------------------------------------------------------------

private:
	Matrix()
	{
		ptr = new double[width * height];
		SetAll(0);
	}

private:
	double* ptr;
};