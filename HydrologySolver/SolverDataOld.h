#pragma once

#include "MultigridMesh.h"
#include <iostream>

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
		B.SetAll(0);
		W.SetAll(0);
		R.SetAll(0);
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

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------

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
				finest.J_pp.Set(i, j, 0);
				finest.J_pm.Set(i, j, 0);
				finest.J_mp.Set(i, j, 0);
				finest.J_mm.Set(i, j, 0);
			}
		}
		finest.UnpackJacobian_From_JVecs();
	}

//-----------------------------------------------------------------------------

	void SolveWithMultigrid(double maxNorm, double minValue, double maxValue, int verbosity, std::filesystem::path outPath)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			finest.F.Scale(-1);
			CalculateJ();
			finest.PrepareMatrixInefficient(outPath);
			while (true)
			{
				finest.VCycleInefficient(maxNorm, 0, 2, outPath);
				std::cout << "V";
				finest.ComputeResidualInefficient();
				double residualNorm = finest.E.Norm();
				if(residualNorm < maxNorm)
					break;
			}

			double correctionNorm = finest.U.Norm();
			if (verbosity > 0)
				std::cout << ", norm: " << correctionNorm << std::endl;
			if (correctionNorm < maxNorm)
				break;
			W.AddCorrectionToHaloedData(finest.U, 1, minValue, maxValue);
			itr++;
		}
	}

//-----------------------------------------------------------------------------

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
			int itrReq = finest.SolveWithGaussSeidelInefficient(maxNorm, 100, verbosity - 1);
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

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------

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


private:
	DataPointer<w + 2, h + 2> B; //Bed Height and boundary condition vectors
	DataPointer<w + 2, h + 2> W; //Water Height and boundary condition vectors
	DataPointer<w, h> R; //Rainfall Values

	MultigridMesh<e> finest;
};

