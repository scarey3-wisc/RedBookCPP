#pragma once
#include "DataPointer.h"
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "MatrixOps.h"

// We assume that all solver data is a square based on a power of 2
// So the dimension will be 2^e + 1

namespace StencilDirections
{
	struct Offset {
		int dx;
		int dy;
	};
	enum Dir {
		cc,
		xp,
		xm,
		yp,
		ym,
		pp,
		pm,
		mp,
		mm
	};
	struct DirList {
		std::array<Dir, 2> dirs;
		uint8_t n;
		const Dir* begin() const { return dirs.data(); }
		const Dir* end() const { return dirs.data() + n; }
	};
	static constexpr std::array<Dir, 9> AllDirs = { cc,xp,xm,yp,ym,pp,pm,mp,mm };
	static constexpr std::array<Dir, 4> Diagonals = { pp,pm,mp,mm };
	static constexpr std::array<Dir, 4> Axes = { xp,xm,yp,ym };
	static constexpr Dir OPPOSITES[9] = {
		cc, //cc
		xm, //xp
		xp, //xm
		ym, //yp
		yp, //ym
		mm, //pp
		mp, //pm
		pm, //mp
		pp  //mm
	};
	static constexpr Dir Opposite(Dir d)
	{
		return OPPOSITES[d];
	}
	static constexpr Offset OFFSETS[9] = {
		{0, 0},   //cc
		{1, 0},   //xp
		{-1, 0},  //xm
		{0, 1},   //yp
		{0, -1},  //ym
		{1, 1},   //pp
		{1, -1},  //pm
		{-1, 1},  //mp
		{-1, -1}  //mm
	};
	static constexpr Dir DIR_FROM_OFF[3][3] = {
	{ Dir::mm, Dir::xm,  Dir::mp },  // dx = -1
	{ Dir::ym, Dir::cc,  Dir::yp  },  // dx =  0
	{ Dir::pm, Dir::xp,  Dir::pp }   // dx = +1
	};
	static constexpr Dir SumDir(Dir a, Dir b)
	{
		Offset oa = OFFSETS[a];
		Offset ob = OFFSETS[b];
		int dx = oa.dx + ob.dx;
		int dy = oa.dy + ob.dy;
		if (dx > 1) dx = 1; else if (dx < 1) dx = -1;
		if (dy > 1) dy = 1; else if (dy < 1) dy = -1;
		return DIR_FROM_OFF[dx + 1][dy + 1];
	}
	static constexpr Offset Del(Dir d)
	{
		return OFFSETS[d];
	}
	static constexpr DirList ORTHOGONAL[9] = {
		{{}, 0},
		{{yp, ym}, 2},
		{{yp, ym}, 2},
		{{xp, xm}, 2},
		{{xp, xm}, 2},
		{{mp, pm}, 2},
		{{mm, pp}, 2},
		{{mm, pp}, 2},
		{{mp, pm}, 2}
	};
	static constexpr DirList Orthogonal(Dir d)
	{
		return ORTHOGONAL[d];
	}
	static constexpr DirList AXIS_ALIGNED[9] = {
		{{}, 0},
		{{xp}, 1},
		{{xm}, 1},
		{{yp}, 1},
		{{ym}, 1},
		{{xp, yp}, 2},
		{{xp, ym}, 2},
		{{xm, yp}, 2},
		{{xm, ym}, 2}
	};
	static constexpr DirList OnAxis(Dir d)
	{
		return AXIS_ALIGNED[d];
	}
	static constexpr DirList DIAG_ALIGNED[9] = {
		{{}, 0},
		{{pp, pm}, 2},
		{{mp, mm}, 2},
		{{pp, mp}, 2},
		{{pm, mm}, 2},
		{{pp}, 1},
		{{pm}, 1},
		{{mp}, 1},
		{{mm}, 1}
	};
	static constexpr DirList OnDiag(Dir d)
	{
		return DIAG_ALIGNED[d];
	}
}

template <int e>
struct MultigridMeshBase
{
	using Dir = StencilDirections::Dir;
	using Offset = StencilDirections::Offset;
	using DirList = StencilDirections::DirList;

	virtual bool VCycle(double maxNorm, int itrDown, int itrUp) = 0;
	virtual bool VCycleInefficient(double maxNorm, int itrDown, int itrUp, std::filesystem::path outPath) = 0;
	virtual void PrepareMatrix(std::filesystem::path outPath) = 0;
	virtual void PrepareMatrixInefficient(std::filesystem::path outPath) = 0;

	static constexpr int w = (1 << e) + 1;
	static constexpr int h = (1 << e) + 1;
	static constexpr int n = w * h;
	static constexpr int wCoarse = (1 << (e - 1)) + 1;
	static constexpr int hCoarse = (1 << (e - 1)) + 1;
	static constexpr int nCoarse = wCoarse * hCoarse;

	constexpr DataPointer<w, h>& JacobianInDir(Dir d)
	{
		switch (d)
		{
		case Dir::cc:
			return J_cc;
		case Dir::xp:
			return J_xp;
		case Dir::xm:
			return J_xm;
		case Dir::yp:
			return J_yp;
		case Dir::ym:
			return J_ym;
		case Dir::pp:
			return J_pp;
		case Dir::pm:
			return J_pm;
		case Dir::mp:
			return J_mp;
		case Dir::mm:
			return J_mm;
		default:
			return J_cc; //should never happen
		}
	}

//-----------------------------------------------------------------------------

	MultigridMeshBase() :
		F(DataPointer<w, h>::AllocateNewPointer()),
		U(DataPointer<w, h>::AllocateNewPointer()),
		E(DataPointer<w, h>::AllocateNewPointer()),
		J_cc(DataPointer<w, h>::AllocateNewPointer()),
		J_xm(DataPointer<w, h>::AllocateNewPointer()),
		J_xp(DataPointer<w, h>::AllocateNewPointer()),
		J_ym(DataPointer<w, h>::AllocateNewPointer()),
		J_yp(DataPointer<w, h>::AllocateNewPointer()),
		J_mm(DataPointer<w, h>::AllocateNewPointer()),
		J_mp(DataPointer<w, h>::AllocateNewPointer()),
		J_pm(DataPointer<w, h>::AllocateNewPointer()),
		J_pp(DataPointer<w, h>::AllocateNewPointer()),
		Jacobian(Matrix<n, n>::AllocateNewPointer()),
		Interpolate(Matrix<nCoarse, n>::AllocateNewPointer()),
		Restrict(Matrix<n, nCoarse>::AllocateNewPointer())
	{
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < h; x++)
			{
				if (x % 2 == 0 && y % 2 == 0)
				{
					Interpolate.Set(x / 2 + y / 2 * wCoarse, x + y * w, 1.0);
				}
				else if (x % 2 == 0 && y % 2 == 1)
				{
					Interpolate.Set(x / 2 + y / 2 * wCoarse, x + y * w, 0.5);
					Interpolate.Set(x / 2 + (y / 2 + 1) * wCoarse, x + y * w, 0.5);
				}
				else if (x % 2 == 1 && y % 2 == 0)
				{
					Interpolate.Set(x / 2 + y / 2 * wCoarse, x + y * w, 0.5);
					Interpolate.Set(x / 2 + 1 + y / 2 * wCoarse, x + y * w, 0.5);
				}
				else
				{
					Interpolate.Set(x / 2 + y / 2 * wCoarse, x + y * w, 0.25);
					Interpolate.Set(x / 2 + 1 + y / 2 * wCoarse, x + y * w, 0.25);
					Interpolate.Set(x / 2 + (y / 2 + 1) * wCoarse, x + y * w, 0.25);
					Interpolate.Set(x / 2 + 1 + (y / 2 + 1) * wCoarse, x + y * w, 0.25);
				}
			}
		}

		MatrixOps::Transpose(Interpolate, Restrict);
		Restrict.Scale(0.25);
	}

	~MultigridMeshBase()
	{
		F.Deallocate();
		U.Deallocate();
		E.Deallocate();
		J_cc.Deallocate();
		J_xm.Deallocate();
		J_xp.Deallocate();
		J_ym.Deallocate();
		J_yp.Deallocate();
		J_mm.Deallocate();
		J_mp.Deallocate();
		J_pm.Deallocate();
		J_pp.Deallocate();
		Jacobian.Deallocate();
		Interpolate.Deallocate();
		Restrict.Deallocate();
	}

	// Disable copying:
	MultigridMeshBase(const MultigridMeshBase&) = delete;
	MultigridMeshBase& operator=(const MultigridMeshBase&) = delete;

//-----------------------------------------------------------------------------

	void
	VisualizeJacobian(std::filesystem::path outPath)
	{
		std::string myString = "Jacobian_" + std::to_string(e) + ".ppm";
		outPath.replace_filename(myString);
		std::ofstream out(outPath.c_str());
		
		out << "P3\n" << w * h << " " << w * h << "\n255\n";
		double min = std::numeric_limits<double>::max();
		double max = std::numeric_limits<double>::min();
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				for (Dir d : StencilDirections::AllDirs)
				{
					double val = JacobianInDir(d).Get(x, y);
					if (val < min)
						min = val;
					if (val > max)
						max = val;
				}
			}
		}
		std::cout << "Jacobian range for e = " << e << ": " << min << " - " << max << std::endl;
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				for (int yy = 0; yy < h; yy++)
				{
					for (int xx = 0; xx < w; xx++)
					{
						double val = 0;
						if (xx == x && yy == y)
							val = J_cc.Get(x, y);
						else if (xx == x + 1 && yy == y)
							val = J_xp.Get(x, y);
						else if (xx == x - 1 && yy == y)
							val = J_xm.Get(x, y);
						else if (xx == x && yy == y + 1)
							val = J_yp.Get(x, y);
						else if (xx == x && yy == y - 1)
							val = J_ym.Get(x, y);
						else if (xx == x + 1 && yy == y + 1)
							val = J_pp.Get(x, y);
						else if (xx == x + 1 && yy == y - 1)
							val = J_pm.Get(x, y);
						else if (xx == x - 1 && yy == y + 1)
							val = J_mp.Get(x, y);
						else if (xx == x - 1 && yy == y - 1)
							val = J_mm.Get(x, y);

						int r = 0;
						int g = 0;
						int b = 0;
						if(val < 0)
							r = val / min * 255;
						else if (val > 0)
							g = val / max * 255;
						out << r << " " << g << " " << b << " ";
					}
				}
				out << "\n";
			}
		}
		out.close();
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationInefficient()
	{
		// Compute an exact solution for U(i, j)
		auto SmoothIndex = [&](int i, int j)
			{
				int index = i + j * w;
				double rhs = F.Get(i, j);
				for (int otherIndex = 0; otherIndex < w * h; otherIndex++)
				{
					if (otherIndex == index)
						continue;
					int oI = otherIndex % w;
					int oJ = otherIndex / w;
					double val = U.Get(oI, oJ);
					val *= Jacobian.Get(otherIndex, index);
					rhs -= val;
				}

				double div = Jacobian.Get(index, index);
				double result = rhs / div;
				if (!std::isnormal(result) && result != 0)
				{
					result = 0;
				}
				U.Set(i, j, result);
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

//-----------------------------------------------------------------------------

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

				if (i > 0 && j > 0)
					rhs -= J_mm.Get(i, j) * U.Get(i - 1, j - 1);
				if (i > 0 && j < h - 1)
					rhs -= J_mp.Get(i, j) * U.Get(i - 1, j + 1);
				if (i < w - 1 && j > 0)
					rhs -= J_pm.Get(i, j) * U.Get(i + 1, j - 1);
				if (i < w - 1 && j < h - 1)
					rhs -= J_pp.Get(i, j) * U.Get(i + 1, j + 1);

				double div = J_cc.Get(i, j);
				double result = rhs / div;
				if (!std::isnormal(result) && result != 0)
				{
					result = 0;
				}
				U.Set(i, j, result);
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

//-----------------------------------------------------------------------------

	int SolveWithGaussSeidelInefficient(double maxNorm, int maxItr, int verbosityLevel)
	{
		U.SetAll(0);
		int itr = 0;
		double prevNorm;
		double initialNorm;
		while (true)
		{
			if (verbosityLevel > 0)
				std::cout << "Gauss Seidel Itr: " << itr;
			ComputeResidualInefficient();
			double norm = E.Norm();
			if (itr == 0)
			{
				prevNorm = norm;
				initialNorm = norm;
			}
			else if (prevNorm < norm && norm > 10 * initialNorm)
			{
				if (verbosityLevel > 0)
					std::cout << " (Diverging, stopping)" << std::endl;
				return -1;
			}
			else
				prevNorm = norm;
			if (verbosityLevel > 0)
				std::cout << ", r: " << norm << std::endl;
			if (norm < maxNorm)
				break;
			GaussSeidelIterationInefficient();
			itr++;
			if (itr > maxItr)
				break;
		}
		return itr;
	}

//-----------------------------------------------------------------------------

	int SolveWithGaussSeidel(double maxNorm, int maxItr, int verbosityLevel)
	{
		U.SetAll(0);
		int itr = 0;
		while (true)
		{
			if (verbosityLevel > 0)
				std::cout << "Gauss Seidel Itr: " << itr;
			ComputeResidual();
			double norm = E.Norm();
			if (verbosityLevel > 0)
				std::cout << ", r: " << norm << std::endl;
			if (norm < maxNorm)
				break;
			GaussSeidelIteration();
			itr++;
			if (itr > maxItr)
				break;
		}
		return itr;
	}

//-----------------------------------------------------------------------------

	void Restrict_Residual_To(MultigridMeshBase<e - 1>& coarse)
	{
		MatrixOps::Multiply(Restrict, E, coarse.F);
	}

//-----------------------------------------------------------------------------

	void Restrict_RHS_From_Residual(MultigridMeshBase<e + 1>& fine)
	{
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				if (i == 0 && j == 0)
					F.Set(i, j, fine.E.Get(2 * i, 2 * j));
				else if (i == 0 && j == h - 1)
					F.Set(i, j, fine.E.Get(2 * i, 2 * j));
				else if (i == w - 1 && j == 0)
					F.Set(i, j, fine.E.Get(2 * i, 2 * j));
				else if (i == w - 1 && j == h - 1)
					F.Set(i, j, fine.E.Get(2 * i, 2 * j));
				else if (i == 0)
				{
					double val = 0;
					val += fine.E.Get(2 * i, 2 * j) / 2;
					val += fine.E.Get(2 * i, 2 * j - 1) / 4;
					val += fine.E.Get(2 * i, 2 * j + 1) / 4;
					F.Set(i, j, val);
				}
				else if (i == w - 1)
				{
					double val = 0;
					val += fine.E.Get(2 * i, 2 * j) / 2;
					val += fine.E.Get(2 * i, 2 * j - 1) / 4;
					val += fine.E.Get(2 * i, 2 * j + 1) / 4;
					F.Set(i, j, val);
				}
				else if (j == 0)
				{
					double val = 0;
					val += fine.E.Get(2 * i, 2 * j) / 2;
					val += fine.E.Get(2 * i + 1, 2 * j) / 4;
					val += fine.E.Get(2 * i - 1, 2 * j) / 4;
					F.Set(i, j, val);
				}
				else if (j == h - 1)
				{
					double val = 0;
					val += fine.E.Get(2 * i, 2 * j) / 2;
					val += fine.E.Get(2 * i + 1, 2 * j) / 4;
					val += fine.E.Get(2 * i - 1, 2 * j) / 4;
					F.Set(i, j, val);
				}
				else
				{
					double val = 0;
					val += fine.E.Get(2 * i, 2 * j) / 4;
					val += fine.E.Get(2 * i + 1, 2 * j) / 8;
					val += fine.E.Get(2 * i - 1, 2 * j) / 8;
					val += fine.E.Get(2 * i, 2 * j + 1) / 8;
					val += fine.E.Get(2 * i, 2 * j - 1) / 8;
					val += fine.E.Get(2 * i + 1, 2 * j + 1) / 16;
					val += fine.E.Get(2 * i - 1, 2 * j + 1) / 16;
					val += fine.E.Get(2 * i + 1, 2 * j - 1) / 16;
					val += fine.E.Get(2 * i - 1, 2 * j - 1) / 16;
					F.Set(i, j, val);
				}
			}
		}
	}

//-----------------------------------------------------------------------------

	void Interpolate_Correction_From(MultigridMeshBase<e - 1>& coarse)
	{
		MatrixOps::Multiply(Interpolate, coarse.U, E);
		U.Add(E, 1.0);
	}

//-----------------------------------------------------------------------------

	void Interpolate_Correction_To(MultigridMeshBase<e + 1>& fine)
	{
		for (int j = 0; j < fine.h; j++)
		{
			for (int i = 0; i < fine.w; i++)
			{
				if (i % 2 == 0 && j % 2 == 0)
				{
					fine.U.Add(i, j, U.Get(i / 2, j / 2));
				}
				else if (i % 2 == 0)
				{
					double val = 0;
					val += U.Get(i / 2, j / 2) / 2;
					val += U.Get(i / 2, j / 2 + 1) / 2;
					fine.U.Add(i, j, val);
				}
				else if (j % 2 == 0)
				{
					double val = 0;
					val += U.Get(i / 2, j / 2) / 2;
					val += U.Get(i / 2 + 1, j / 2) / 2;
					fine.U.Add(i, j, val);
				}
				else
				{
					double val = 0;
					val += U.Get(i / 2, j / 2) / 4;
					val += U.Get(i / 2 + 1, j / 2) / 4;
					val += U.Get(i / 2, j / 2 + 1) / 4;
					val += U.Get(i / 2 + 1, j / 2 + 1) / 4;
					fine.U.Add(i, j, val);
				}
			}
		}
	}

//-----------------------------------------------------------------------------

	void UnpackJacobian_From_JVecs()
	{
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int index = x + y * w;
				Jacobian.Set(index, index, J_cc.Get(x, y));
				if (x < w - 1)
					Jacobian.Set(index + 1, index, J_xp.Get(x, y));
				if (x > 0)
					Jacobian.Set(index - 1, index, J_xm.Get(x, y));
				if (y < h - 1)
					Jacobian.Set(index + w, index, J_yp.Get(x, y));
				if (y > 0)
					Jacobian.Set(index - w, index, J_ym.Get(x, y));
				if (x < w - 1 && y < h - 1)
					Jacobian.Set(index + 1 + w, index, J_pp.Get(x, y));
				if (x < w - 1 && y > 0)
					Jacobian.Set(index + 1 - w, index, J_pm.Get(x, y));
				if (x > 0 && y < h - 1)
					Jacobian.Set(index - 1 + w, index, J_mp.Get(x, y));
				if (x > 0 && y > 0)
					Jacobian.Set(index - 1 - w, index, J_mm.Get(x, y));
			}
		}
	}

//-----------------------------------------------------------------------------

	void Pack_Jacobian_In_JVecs()
	{
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int index = x + y * w;
				J_cc.Set(x, y, Jacobian.Get(index, index));
				if (x < w - 1)
					J_xp.Set(x, y, Jacobian.Get(index + 1, index));
				if (x > 0)
					J_xm.Set(x, y, Jacobian.Get(index - 1, index));
				if (y < h - 1)
					J_yp.Set(x, y, Jacobian.Get(index + w, index));
				if (y > 0)
					J_ym.Set(x, y, Jacobian.Get(index - w, index));
				if (x < w - 1 && y < h - 1)
					J_pp.Set(x, y, Jacobian.Get(index + 1 + w, index));
				if (x < w - 1 && y > 0)
					J_pm.Set(x, y, Jacobian.Get(index + 1 - w, index));
				if (x > 0 && y < h - 1)
					J_mp.Set(x, y, Jacobian.Get(index - 1 + w, index));
				if (x > 0 && y > 0)
					J_mm.Set(x, y, Jacobian.Get(index - 1 - w, index));
			}
		}
	}

//-----------------------------------------------------------------------------

	void Restrict_Jacobian_To(MultigridMeshBase<e - 1>& coarse)
	{
		Matrix<n, nCoarse> temp = Matrix<n, nCoarse>::AllocateNewPointer();
		MatrixOps::Multiply(Restrict, Jacobian, temp);
		MatrixOps::Multiply(temp, Interpolate, coarse.Jacobian);
		temp.Deallocate();
	}

//-----------------------------------------------------------------------------

	void Restrict_Jacobian_From(MultigridMeshBase<e + 1>& fine)
	{
		J_xp.SetAll(0);
		J_xm.SetAll(0);
		J_yp.SetAll(0);
		J_ym.SetAll(0);
		J_pp.SetAll(0);
		J_pm.SetAll(0);
		J_mp.SetAll(0);
		J_mm.SetAll(0);
		J_cc.SetAll(0);
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				if (i == 0 && j == 0)
					Restrict_Jacobian_On_Corner_From<Dir::pp>(fine, i, j);
				else if (i == 0 && j == h - 1)
					Restrict_Jacobian_On_Corner_From<Dir::pm>(fine, i, j);
				else if (i == w - 1 && j == 0)
					Restrict_Jacobian_On_Corner_From<Dir::mp>(fine, i, j);
				else if (i == w - 1 && j == h - 1)
					Restrict_Jacobian_On_Corner_From<Dir::mm>(fine, i, j);
				else if (i == 0)
					Restrict_Jacobian_On_Side_From<Dir::xp>(fine, i, j);
				else if (i == w - 1)
					Restrict_Jacobian_On_Side_From<Dir::xm>(fine, i, j);
				else if (j == 0)
					Restrict_Jacobian_On_Side_From<Dir::yp>(fine, i, j);
				else if (j == h - 1)
					Restrict_Jacobian_On_Side_From<Dir::ym>(fine, i, j);
				else
					Restrict_Jacobian_On_Interior_From(fine, i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	template<Dir inward>
	inline void Restrict_Jacobian_On_Corner_From(MultigridMeshBase<e + 1>& fine, int i, int j)
	{
		double selfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i, 2 * j);
		JacobianInDir(Dir::cc).Add(i, j, selfTerm);

		for (Dir a : OnAxis(inward))
		{
			double adjTerm = fine.JacobianInDir(a).Get(2 * i, 2 * j) / 2;
			JacobianInDir(Dir::cc).Add(i, j, adjTerm);
			JacobianInDir(a).Add(i, j, adjTerm);
		}

		double diagTerm = fine.JacobianInDir(inward).Get(2 * i, 2 * j) / 4;
		JacobianInDir(Dir::cc).Add(i, j, diagTerm);
		JacobianInDir(inward).Add(i, j, diagTerm);
		for (Dir a : OnAxis(inward))
		{
			JacobianInDir(a).Add(i, j, diagTerm);
		}
	}

//-----------------------------------------------------------------------------

	template<Dir inward>
	void Restrict_Jacobian_On_Side_From(MultigridMeshBase<e + 1>& fine, int i, int j)
	{
		double selfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i, 2 * j) / 2;
		JacobianInDir(Dir::cc).Add(i, j, selfTerm);

		double inTerm = fine.JacobianInDir(inward).Get(2 * i, 2 * j) / 4;
		JacobianInDir(Dir::cc).Add(i, j, inTerm);
		JacobianInDir(inward).Add(i, j, inTerm);

		for (Dir o : Orthogonal(inward))
		{
			double orthoTerm = fine.JacobianInDir(o).Get(2 * i, 2 * j) / 4;
			JacobianInDir(Dir::cc).Add(i, j, orthoTerm);
			JacobianInDir(o).Add(i, j, orthoTerm);
		}

		for (Dir d : OnDiag(inward))
		{
			double diagTerm = fine.JacobianInDir(d).Get(2 * i, 2 * j) / 8;
			JacobianInDir(Dir::cc).Add(i, j, diagTerm);
			JacobianInDir(d).Add(i, j, diagTerm);
			for (Dir o : Orthogonal(d))
			{
				JacobianInDir(o).Add(i, j, diagTerm);
			}
		}

		for (Dir o : Orthogonal(inward))
		{
			Offset del = Del(o);
			double oSelfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i + del.dx, 2 * j + del.dy) / 8;
			JacobianInDir(Dir::cc).Add(i, j, oSelfTerm);
			JacobianInDir(o).Add(i, j, oSelfTerm);

			double oInTerm = fine.JacobianInDir(inward).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
			JacobianInDir(Dir::cc).Add(i, j, oInTerm);
			JacobianInDir(o).Add(i, j, oInTerm);
			JacobianInDir(inward).Add(i, j, oInTerm);
			JacobianInDir(SumDir(o, inward)).Add(i, j, oInTerm);

			double oForwardTerm = fine.JacobianInDir(o).Get(2 * i + del.dx, 2 * j + del.dy) / 4;
			JacobianInDir(o).Add(i, j, oForwardTerm);

			double oBackwardTerm = fine.JacobianInDir(Opposite(o)).Get(2 * i + del.dx, 2 * j + del.dy) / 4;
			JacobianInDir(Dir::cc).Add(i, j, oBackwardTerm);

			double oForwardInTerm = fine.JacobianInDir(SumDir(o, inward)).Get(2 * i + del.dx, 2 * j + del.dy) / 8;
			JacobianInDir(o).Add(i, j, oForwardInTerm);
			JacobianInDir(SumDir(o, inward)).Add(i, j, oForwardInTerm);

			double oBackwardInTerm = fine.JacobianInDir(SumDir(Opposite(o), inward)).Get(2 * i + del.dx, 2 * j + del.dy) / 8;
			JacobianInDir(Dir::cc).Add(i, j, oBackwardInTerm);
			JacobianInDir(inward).Add(i, j, oBackwardInTerm);
		}
	}

//-----------------------------------------------------------------------------

	void Restrict_Jacobian_On_Interior_From(MultigridMeshBase<e + 1>& fine, int i, int j)
	{
		double selfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i, 2 * j) / 4;
		JacobianInDir(Dir::cc).Add(i, j, selfTerm);

		for (Dir a : StencilDirections::Axes)
		{
			double adjTerm = fine.JacobianInDir(a).Get(2 * i, 2 * j) / 8;
			JacobianInDir(Dir::cc).Add(i, j, adjTerm);
			JacobianInDir(a).Add(i, j, adjTerm);
		}

		for (Dir d : StencilDirections::Diagonals)
		{
			double diagTerm = fine.JacobianInDir(d).Get(2 * i, 2 * j) / 16;
			JacobianInDir(Dir::cc).Add(i, j, diagTerm);
			JacobianInDir(d).Add(i, j, diagTerm);
			for (Dir a : OnAxis(d))
			{
				JacobianInDir(a).Add(i, j, diagTerm);
			}
		}

		for (Dir a : StencilDirections::Axes)
		{
			Offset del = Del(a);
			double adjSelfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
			JacobianInDir(Dir::cc).Add(i, j, adjSelfTerm);
			JacobianInDir(a).Add(i, j, adjSelfTerm);

			double adjAdjTerm = fine.JacobianInDir(a).Get(2 * i + del.dx, 2 * j + del.dy) / 8;
			JacobianInDir(a).Add(i, j, adjAdjTerm);

			double adjRevTerm = fine.JacobianInDir(Opposite(a)).Get(2 * i + del.dx, 2 * j + del.dy) / 8;
			JacobianInDir(Dir::cc).Add(i, j, adjRevTerm);

			for (Dir adjO : Orthogonal(a))
			{
				double adjOrtTerm = fine.JacobianInDir(adjO).Get(2 * i + del.dx, 2 * j + del.dy) / 32;

				JacobianInDir(Dir::cc).Add(i, j, adjOrtTerm);
				JacobianInDir(adjO).Add(i, j, adjOrtTerm);
				JacobianInDir(a).Add(i, j, adjOrtTerm);
				JacobianInDir(SumDir(a, adjO)).Add(i, j, adjOrtTerm);
			}
			for (Dir adjAdjDiag : OnDiag(a))
			{
				double adjAdjDiagTerm = fine.JacobianInDir(adjAdjDiag).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
				JacobianInDir(a).Add(i, j, adjAdjDiagTerm);
				JacobianInDir(adjAdjDiag).Add(i, j, adjAdjDiagTerm);
			}
			for (Dir adjO : Orthogonal(a))
			{
				Dir adjOppDiag = SumDir(Opposite(a), adjO);
				double adjOppDiagTerm = fine.JacobianInDir(adjOppDiag).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
				JacobianInDir(Dir::cc).Add(i, j, adjOppDiagTerm);
				JacobianInDir(adjO).Add(i, j, adjOppDiagTerm);
			}
		}

		for (Dir d : StencilDirections::Diagonals)
		{
			Offset del = Del(d);
			double diagSelfTerm = fine.JacobianInDir(Dir::cc).Get(2 * i + del.dx, 2 * j + del.dy) / 64;
			JacobianInDir(Dir::cc).Add(i, j, diagSelfTerm);
			JacobianInDir(d).Add(i, j, diagSelfTerm);
			for (Dir a : OnAxis(d))
			{
				JacobianInDir(a).Add(i, j, diagSelfTerm);
			}

			double diagDiagTerm = fine.JacobianInDir(d).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
			JacobianInDir(d).Add(i, j, diagDiagTerm);

			double diagRevTerm = fine.JacobianInDir(Opposite(d)).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
			JacobianInDir(Dir::cc).Add(i, j, diagRevTerm);

			for (Dir o : Orthogonal(d))
			{
				double diagOrtTerm = fine.JacobianInDir(o).Get(2 * i + del.dx, 2 * j + del.dy) / 16;
				JacobianInDir(SumDir(o, d)).Add(i, j, diagOrtTerm);
			}

			for (Dir a : OnAxis(d))
			{
				double diagAxisForwardTerm = fine.JacobianInDir(a).Get(2 * i + del.dx, 2 * j + del.dy) / 32;
				JacobianInDir(a).Add(i, j, diagAxisForwardTerm);
				JacobianInDir(d).Add(i, j, diagAxisForwardTerm);
			}

			for (Dir o : OnAxis(Opposite(d)))
			{
				double diagAxisBackwardTerm = fine.JacobianInDir(o).Get(2 * i + del.dx, 2 * j + del.dy) / 32;
				JacobianInDir(Dir::cc).Add(i, j, diagAxisBackwardTerm);
				JacobianInDir(SumDir(o, d)).Add(i, j, diagAxisBackwardTerm);
			}
		}
	}

//-----------------------------------------------------------------------------

	void ComputeResidualInefficient()
	{
		MatrixOps::Multiply(Jacobian, U, E);
		E.Add(F, -1);
		E.Scale(-1);
	}

//-----------------------------------------------------------------------------

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
				if (i > 0 && j > 0)
					result += J_mm.Get(i, j) * U.Get(i - 1, j - 1);
				if (i > 0 && j < h - 1)
					result += J_mp.Get(i, j) * U.Get(i - 1, j + 1);
				if (i < w - 1 && j > 0)
					result += J_pm.Get(i, j) * U.Get(i + 1, j - 1);
				if (i < w - 1 && j < h - 1)
					result += J_pp.Get(i, j) * U.Get(j + 1, j + 1);

				result -= F.Get(i, j);
				result *= -1;
				E.Set(i, j, result);
			}
		}
	}

//-----------------------------------------------------------------------------

	DataPointer<w, h> F; //buffer for evaluations of F_ij(W)
	DataPointer<w, h> U; //buffer for correction that we calculate
	DataPointer<w, h> E; //buffer for residual after evaluating AU

	/*
	* The Jacobian of F is a pentadiagonal matrix with each diagonal
	* corresponding to one spoke of the 5 point stencil.
	* We store each diagonal (each spoke of the stencil) in a separate
	* DataPointer for maximum efficiency in sparsity
	*
	* As we use multigrid, we might also populate the diagonals of
	* coarser meshes, and so our 5 point stencil representation may
	* become a 9 point stencil on coarser meshes. For that purpose,
	* we allucate 4 more diagonals:
	*/
	DataPointer<w, h> J_cc; //Center
	DataPointer<w, h> J_xp; //i+1, j
	DataPointer<w, h> J_xm; //i-1, j
	DataPointer<w, h> J_yp; //i, j+1
	DataPointer<w, h> J_ym; //i, j-1
	DataPointer<w, h> J_mm; //i-1, j-1
	DataPointer<w, h> J_mp; //i-1, j+1
	DataPointer<w, h> J_pm; //i+1, j-1
	DataPointer<w, h> J_pp; //i+1, j+1

	Matrix<n, n> Jacobian;
	Matrix<nCoarse, n> Interpolate;
	Matrix<n, nCoarse> Restrict;

};


template<int e>
struct MultigridMesh : public MultigridMeshBase<e>{
	static_assert(e > 4);

	MultigridMesh<e - 1> child;
	bool VCycle(double maxNorm, int itrDown, int itrUp)
	{
		for (int i = 0; i < itrDown; i++)
			this->GaussSeidelIteration();
		this->ComputeResidual();
		child.Restrict_RHS_From_Residual(*this);
		child.U.SetAll(0);
		child.VCycle(maxNorm, itrDown, itrUp);
		child.Interpolate_Correction_To(*this);
		for (int i = 0; i < itrUp; i++)
			this->GaussSeidelIteration();
		return true;
	}
	bool VCycleInefficient(double maxNorm, int itrDown, int itrUp, std::filesystem::path outPath)
	{
		MatrixOps::Visualize(this->F, outPath, "FFull_" + std::to_string(e));
		for (int i = 0; i < itrDown; i++)
			this->GaussSeidelIterationInefficient();
		this->ComputeResidualInefficient();
		this->Restrict_Residual_To(child);
		child.U.SetAll(0);
		bool okay = child.VCycleInefficient(maxNorm, itrDown, itrUp, outPath);
		if (okay)
		{
			this->Interpolate_Correction_From(child);
			for (int i = 0; i < itrUp; i++)
				this->GaussSeidelIterationInefficient();
			return true;
		}
		else
		{
			return this->SolveWithGaussSeidelInefficient(maxNorm, 10000, 0);
		}

	}
	void PrepareMatrix(std::filesystem::path outPath)
	{
		this->VisualizeJacobian(outPath);
		child.Restrict_Jacobian_From(*this);
		child.PrepareMatrix(outPath);
	}
	void PrepareMatrixInefficient(std::filesystem::path outPath)
	{
		MatrixOps::Visualize(this->Jacobian, outPath, "JacobianFull_" + std::to_string(e));
		this->Restrict_Jacobian_To(child);
		child.PrepareMatrixInefficient(outPath);
	}
};

template <>
struct MultigridMesh<4> : public MultigridMeshBase<4>
{
	bool VCycle(double maxNorm, int itrDown, int itrUp)
	{
		SolveWithGaussSeidel(maxNorm, 100, 0);
		return true;
	}
	bool VCycleInefficient(double maxNorm, int itrDown, int itrUp, std::filesystem::path outPath)
	{
		MatrixOps::Visualize(this->F, outPath, "FFull_" + std::to_string(4));
		int itr = SolveWithGaussSeidelInefficient(maxNorm, 10000, 0);
		if(itr < 0)
			return false;
		else
			return true;
	}
	void PrepareMatrix(std::filesystem::path outPath)
	{
		this->VisualizeJacobian(outPath);
	}
	void PrepareMatrixInefficient(std::filesystem::path outPath)
	{
		MatrixOps::Visualize(this->Jacobian, outPath, "JacobianFull_" + std::to_string(4));
	}
};