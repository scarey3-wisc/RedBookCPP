#pragma once

#include "DataPointer.h"
#include <iostream>

// We assume that all solver data is a square based on a power of 2
// So the dimension will be 2^e + 1

template <int e>
class MultigridMesh
{
public:
	MultigridMesh() :
		F(DataPointer<w, h>::AllocateNewPointer()),
		U(DataPointer<w, h>::AllocateNewPointer()),
		E(DataPointer<w, h>::AllocateNewPointer()),
		J_cc(DataPointer<w, h>::AllocateNewPointer()),
		J_xm(DataPointer<w, h>::AllocateNewPointer()),
		J_xp(DataPointer<w, h>::AllocateNewPointer()),
		J_ym(DataPointer<w, h>::AllocateNewPointer()),
		J_yp(DataPointer<w, h>::AllocateNewPointer())
	{}

	~MultigridMesh()
	{
		F.Deallocate();
		U.Deallocate();
		J_cc.Deallocate();
		J_xm.Deallocate();
		J_xp.Deallocate();
		J_ym.Deallocate();
		J_yp.Deallocate();
	}

	// Disable copying:
	MultigridMesh(const MultigridMesh&) = delete;
	MultigridMesh& operator=(const MultigridMesh&) = delete;

	void RestrictRHSTo(MultigridMesh<e - 1>& coarser);
	void InterpolateUTo(MultigridMesh<e + 1>& finer);

	void GaussSeidelIteration()
	{
		// Compute an exact solution for U(i, j)
		auto SmoothIndex = [&](int i, int j)
			{
				double rhs = F.Get(i, j);
				if (i > 0)
					rhs -= J_xm.Get(i, j) * U.Get(i - 1, j);
				if (i < w - 1)
					rhs -= J_xp.Get(i, j) * U.Get(i + 1, j);
				if (j > 0)
					rhs -= J_ym.Get(i, j) * U.Get(i, j - 1);
				if (j < h - 1)
					rhs -= J_yp.Get(i, j) * U.Get(i, j + 1);
				double div = J_cc.Get(i, j);
				U.Set(i, j, rhs / div);
			};

		//Red Iteration
		for (int j = 0; j < h; j++)
		{
			for (int i = j % 2; i < w; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
		//Black Iteration
		for (int j = 0; j < h; j++)
		{
			for (int i = 1 - j % 2; i < w; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}
	int SolveWithGaussSeidel(double maxNorm, int verbosityLevel)
	{
		U.Set(0);
		int itr = 0;
		while (true)
		{
			if (verbosityLevel > 0)
				std::cout << "Gauss Seidel Itr: " << itr;
			ComputeResidual();
			double norm = E.Norm();
			if(verbosityLevel > 0)
				std::cout << ", r: " << norm << std::endl;
			if (norm < maxNorm)
				break;
			GaussSeidelIteration();
			itr++;
		}
		return itr;
	}

	static constexpr int w = (1 << e) + 1;
	static constexpr int h = (1 << e) + 1;

private:
	void ComputeResidual()
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				double result = J_cc.Get(i, j) * U.Get(i, j);
				if (i > 0)
					result += J_xm.Get(i, j) * U.Get(i - 1, j);
				if (i < w - 1)
					result += J_xp.Get(i, j) * U.Get(i + 1, j);
				if (j > 0)
					result += J_ym.Get(i, j) * U.Get(i, j - 1);
				if (j < h - 1)
					result += J_yp.Get(i, j) * U.Get(i, j + 1);

				result -= F.Get(i, j);
				E.Set(i, j, result);
			}
		}
	}

public:
	DataPointer<w, h> F; //buffer for evaluations of F_ij(W)
	DataPointer<w, h> U; //buffer for correction that we calculate
	DataPointer<w, h> E; //buffer for residual after evaluating AU

	/*
	* The Jacobian of F is a pentadiagonal matrix with each diagonal
	* corresponding to one spoke of the 5 point stencil.
	* We store each diagonal (each spoke of the stencil) in a separate
	* DataPointer for maximum efficiency in sparsity
	* So J_cc corresponds with the central part of the stencil
	* J_xp corresponds to the i+1, j part of the stencil
	* J_xm corresponds to the i-1, j part of the stencil
	* J_yp corresponds to the i, j+1 part of the stencil
	* J_ym corresponds to the i, j-1 part of the stencil
	*
	* Note that some of these values may be garabage - for instance,
	* the (0, 5) grid location won't have a spoke in the -x direction,
	* but we will have allocated space for it.
	* It'll be the responsibility of the Jacobian Calculator to ensure
	* that we don't waste time calculating these values, and of the
	* Newton Solver to ensure that we don't use them.
	*/
	DataPointer<w, h> J_cc;
	DataPointer<w, h> J_xp;
	DataPointer<w, h> J_xm;
	DataPointer<w, h> J_yp;
	DataPointer<w, h> J_ym;
};

//-----------------------------------------------------------------------------
//*****************************************************************************
//-----------------------------------------------------------------------------

template <int e>
class SolverData
{
public:
	SolverData(double gravity, double drag, double spacing) :
		g(gravity), l(drag), s(spacing), fc(0.5 * gravity / drag), s2(spacing * spacing),
		B(DataPointer<w + 2, h + 2>::AllocateNewPointer()),
		W(DataPointer<w + 2, h + 2>::AllocateNewPointer()),
		R(DataPointer<w, h>::AllocateNewPointer())
	{

	}
	~SolverData()
	{
		B.Deallocate();
		W.Deallocate();
		R.Deallocate();
	}

	// Disable copying:
	SolverData(const SolverData&) = delete;
	SolverData& operator=(const SolverData&) = delete;

	void CalculateF()
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				double del = GetR(i, j);
				del -= Flux(i, j, i + 1, j) / s2;
				del -= Flux(i, j, i - 1, j) / s2;
				del -= Flux(i, j, i, j + 1) / s2;
				del -= Flux(i, j, i, j - 1) / s2;
				finest.F.Set(i, j, del);
			}
		}
	}
	void CalculateJ()
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				double daxp = DFluxDA(i, j, i + 1, j);
				double daxm = DFluxDA(i, j, i - 1, j);
				double dayp = DFluxDA(i, j, i, j + 1);
				double daym = DFluxDA(i, j, i, j - 1);
				finest.J_cc.Set(i, j, -1 * (daxp + daxm + dayp + daym) / s2);
				finest.J_xp.Set(i, j, -1 * DFluxDB(i, j, i + 1, j) / s2);
				finest.J_xm.Set(i, j, -1 * DFluxDB(i, j, i - 1, j) / s2);
				finest.J_yp.Set(i, j, -1 * DFluxDB(i, j, i, j + 1) / s2);
				finest.J_ym.Set(i, j, -1 * DFluxDB(i, j, i, j - 1) / s2);
			}
		}
	}

	void SolveWithGaussSeidel(double maxNorm, double minValue, double maxValue, int verbosity)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			finest.F.Scale(-1);
			CalculateJ();
			int itrReq = finest.SolveWithGaussSeidel(maxNorm, verbosity - 1);
			if (itrReq > 0)
				std::cout << ", requiring " << itrReq << " GS iterations";
			double correctionNorm = finest.U.Norm();
			if (verbosity > 0)
				std::cout << ", norm: " << correctionNorm << std::endl;
			if (correctionNorm < maxNorm)
				break;
			W.AddCorrectionToHaloedData(finest.U, 1, minValue, maxValue);
			itr++;
		}
	}

	inline double GetW(int x, int y)
	{
		return W.Get(x + 1, y + 1);
	}
	inline double GetZ(int x, int y)
	{
		return GetB(x, y) + GetW(x, y);
	}
	inline double GetB(int x, int y)
	{
		return B.Get(x + 1, y + 1);
	}
	inline double GetR(int x, int y)
	{
		return R.Get(x, y);
	}
	inline void SetW(int x, int y, double value)
	{
		W.Set(x + 1, y + 1, value);
	}
	inline void SetB(int x, int y, double value)
	{
		B.Set(x + 1, y + 1, value);
	}
	inline void SetR(int x, int y, double value)
	{
		R.Set(x, y, value);
	}

	static constexpr int w = (1 << e) + 1;
	static constexpr int h = (1 << e) + 1;

	const double g; // gravitational constant
	const double l; // drag constant
	const double s; // spacing constant

	const double fc; // the quantity g/2l, an important flux constant
	const double s2; // spacing squared is often an important quantity
private:

//-----------------------------------------------------------------------------
// Returns a positive value if flux is flowing from (ax, ay) to (bx, by)
	double Flux(int ax, int ay, int bx, int by)
	{
		double zDiff = GetZ(ax, ay) - GetZ(bx, by);
		if (zDiff > 0)
			return fc * zDiff * GetW(ax, ay);
		else if (zDiff < 0)
			return fc * zDiff * GetW(bx, by);
		else
			return 0;
	}

//-----------------------------------------------------------------------------
// Derivative in flux between (ax, ay) to (bx, by) with respect to water in A
	double DFluxDA(int ax, int ay, int bx, int by)
	{
		double zDiff = GetZ(ax, ay) - GetZ(bx, by);
		if (zDiff > 0)
			return fc * zDiff + fc * GetW(ax, ay);
		else if (zDiff < 0)
			return fc * GetW(bx, by);
		else
			return fc * (GetW(ax, ay) + GetW(bx, by)) / 2;
	}

//-----------------------------------------------------------------------------
// Derivative in flux between (ax, ay) to (bx, by) with respect to water in B
	double DFluxDB(int ax, int ay, int bx, int by)
	{
		double zDiff = GetZ(ax, ay) - GetZ(bx, by);
		if (zDiff > 0)
			return -1 * fc * GetW(ax, ay);
		else if (zDiff < 0)
			return -1 * fc * GetW(bx, by) + fc * zDiff;
		else
			return -fc * (GetW(ax, ay) + GetW(bx, by)) / 2;
	}

	DataPointer<w + 2, h + 2> B; //Bed Height and boundary condition vectors
	DataPointer<w + 2, h + 2> W; //Water Height and boundary condition vectors
	DataPointer<w, h> R; //Rainfall Values

	MultigridMesh<e> finest;
};

