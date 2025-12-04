#pragma once

#pragma once

#include "MatrixOps.h"
#include "CSRMatrix.h"
#include "mkl.h"
#include <iostream>
#include <omp.h>

template <int e>
class SolverDataBase
{
public:
	SolverDataBase(double gravity, double drag, double spacing) :
		g(gravity), l(drag), s(spacing), fc(0.5 * gravity / drag), s2(spacing* spacing),
		B(DataPointer<w, h>::AllocateNewPointer()),
		W(DataPointer<w, h>::AllocateNewPointer()),
		R(DataPointer<w, h>::AllocateNewPointer()),
		F(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		U(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		O(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		J_cc(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		J_xp(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		J_xm(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		J_yp(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		J_ym(DataPointer<w - 2, h - 2>::AllocateNewPointer()),
		JS(FiveStencilCSRMatrix<w - 2, h - 2>::AllocateNewPointer())
	{

	}
	~SolverDataBase()
	{
		B.Deallocate();
		W.Deallocate();
		R.Deallocate();
		F.Deallocate();
		U.Deallocate();
		O.Deallocate();
		J_cc.Deallocate();
		J_xp.Deallocate();
		J_xm.Deallocate();
		J_yp.Deallocate();
		J_ym.Deallocate();
		JS.Deallocate();
	}

	// Disable copying:
	SolverDataBase(const SolverDataBase&) = delete;
	SolverDataBase& operator=(const SolverDataBase&) = delete;

	//-----------------------------------------------------------------------------

	void CalculateF()
	{
#pragma omp parallel for
		for (int j = 1; j < h - 1; j++)
		{
			for (int i = 1; i < w - 1; i++)
			{
				double del = GetR(i, j);
				del -= Flux(i, j, i + 1, j) / s2;
				del -= Flux(i, j, i - 1, j) / s2;
				del -= Flux(i, j, i, j + 1) / s2;
				del -= Flux(i, j, i, j - 1) / s2;
				F.Set(i - 1, j - 1, del);
			}
		}
	}

//-----------------------------------------------------------------------------

	void CalculateJ()
	{
#pragma omp parallel for
		for (int j = 1; j < h - 1; j++)
		{
			for (int i = 1; i < w - 1; i++)
			{
				double daxp = DFluxDA(i, j, i + 1, j);
				double daxm = DFluxDA(i, j, i - 1, j);
				double dayp = DFluxDA(i, j, i, j + 1);
				double daym = DFluxDA(i, j, i, j - 1);
				J_cc.Set(i - 1, j - 1, -1 * (daxp + daxm + dayp + daym) / s2);
				J_xp.Set(i - 1, j - 1, -1 * DFluxDB(i, j, i + 1, j) / s2);
				J_xm.Set(i - 1, j - 1, -1 * DFluxDB(i, j, i - 1, j) / s2);
				J_yp.Set(i - 1, j - 1, -1 * DFluxDB(i, j, i, j + 1) / s2);
				J_ym.Set(i - 1, j - 1, -1 * DFluxDB(i, j, i, j - 1) / s2);
			}
		}
	}

//-----------------------------------------------------------------------------

	void SolveSparseSystemWithPARDISO(int matrixSize, double* vals, int* rowOffsets, int* colIndices, double* x, double* f, bool verbose)
	{


		//std::vector<int> perm(262144);
		//for (int i = 0; i < 262144; i++)
		//    perm[i] = i;

		MKL_INT n = matrixSize; // Matrix size
		MKL_INT mtype = 11;        // Real and structurally symmetric matrix
		MKL_INT nrhs = 1;         // Number of right hand sides
		void* pt[64];             // Internal solver memory pointer pt
		// should be "int" when using 32-bit architectures, or "long int"
		// for 64-bit architectures. void* should be OK in both cases
		MKL_INT iparm[64];        // Pardiso control parameters
		MKL_INT maxfct, mnum, phase, error, msglvl;
		MKL_INT i;
		float ddum;               // Scalar dummy (PARDISO needs it)
		MKL_INT idum;             // Integer dummy (PARDISO needs it)

		// Set-up PARDISO control parameters

		for (i = 0; i < 64; i++)
		{
			iparm[i] = 0;
		}
		iparm[0] = 1;         // No solver default
		iparm[1] = 3;         // Parallel nested dissection
		iparm[3] = 0;         // No iterative-direct algorithm
		iparm[4] = 0;         // No user fill-in reducing permutation
		iparm[5] = 0;         // Write solution into x
		iparm[6] = 0;         // Not in use
		iparm[7] = 0;         // Max numbers of iterative refinement steps
		iparm[8] = 0;         // Not in use
		iparm[9] = 13;        // Perturb the pivot elements with 1E-13
		iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS
		iparm[11] = 0;        // Not in use
		iparm[12] = 0;        // Maximum weighted matching algorithm is switched-off
		// Try iparm[12] = 1 in case of inappropriate accuracy
		iparm[13] = 0;        // Output: Number of perturbed pivots
		iparm[14] = 0;        // Not in use
		iparm[15] = 0;        // Not in use
		iparm[16] = 0;        // Not in use
		iparm[17] = -1;       // Output: Number of nonzeros in the factor LU
		iparm[18] = -1;       // Output: Mflops for LU factorization
		iparm[19] = 0;        // Output: Numbers of CG Iterations
		iparm[23] = 1;        // Two-level factorization*/
		iparm[26] = 1;        // Check matrix for errors
		iparm[27] = 0;        // Use double precision
		iparm[34] = 1;        // Use zero-based indexing
		maxfct = 1;           // Maximum number of numerical factorizations.
		mnum = 1;             // Which factorization to use.
		msglvl = verbose? 1 : 0;           // Print statistical information in file
		error = 0;            // Initialize error flag

		// Initialize the internal solver memory pointer. This is only
		// necessary for the FIRST call of the PARDISO solver
		for (i = 0; i < 64; i++)
		{
			pt[i] = 0;
		}

		// Reordering and Symbolic Factorization. This step also allocates
		// all memory that is necessary for the factorization
		phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
			vals, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0)
			throw std::runtime_error("PARDISO error during symbolic factorization");

		if (verbose)
		{
			std::cout << "Reordering completed ... " << std::endl;
			std::cout << "Number of nonzeros in factors = " << iparm[17] << std::endl;
			std::cout << "Number of factorization MFLOPS = " << iparm[18] << std::endl;
		}


		// Numerical factorization
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
			vals, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0)
			throw std::runtime_error("PARDISO error during numerical factorization");

		if (verbose)
		{
			std::cout << "Factorization completed ... " << std::endl;
		}


		// Back substitution and iterative refinement
		phase = 33;
		iparm[7] = 0;         // Max numbers of iterative refinement steps
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
			vals, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, static_cast<void*>(f), x, &error);
		if (error != 0)
			throw std::runtime_error("PARDISO error during solution phase");

		if (verbose)
		{
			std::cout << "Solve completed ... " << std::endl;
		}


		// Termination and release of memory.
		phase = -1;           // Release internal memory
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
			&ddum, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	}


//-----------------------------------------------------------------------------

	void ComputeResidualIntoO()
	{
#pragma omp parallel for
		for (int j = 0; j < h - 2; j++)
		{
			for (int i = 0; i < w - 2; i++)
			{
				double result = J_cc.Get(i, j) * U.Get(i, j);
				if (i > 0)
					result += J_xm.Get(i, j) * U.Get(i - 1, j);
				if (i < w - 3)
					result += J_xp.Get(i, j) * U.Get(i + 1, j);
				if (j > 0)
					result += J_ym.Get(i, j) * U.Get(i, j - 1);
				if (j < h - 3)
					result += J_yp.Get(i, j) * U.Get(i, j + 1);

				result -= F.Get(i, j);
				result *= -1;
				O.Set(i, j, result);
			}
		}
	}

	void SmoothIndex(int i, int j)
	{
		double rhs = F.Get(i, j);
		if (i > 0)
			rhs -= J_xm.Get(i, j) * U.Get(i - 1, j);
		if (i < w - 3)
			rhs -= J_xp.Get(i, j) * U.Get(i + 1, j);
		if (j > 0)
			rhs -= J_ym.Get(i, j) * U.Get(i, j - 1);
		if (j < h - 3)
			rhs -= J_yp.Get(i, j) * U.Get(i, j + 1);


		double div = J_cc.Get(i, j);
		double result = rhs / div;
		if (!std::isnormal(result) && result != 0)
		{
			result = 0;
		}
		U.Set(i, j, result);
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationOnSemiOddPoints()
	{
#pragma omp parallel for
		for (int j = 0; j < h - 2; j+= 2)
		{
			for (int i = 1; i < w - 2; i+= 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		for (int j = 1; j < h - 2; j += 2)
		{
			for (int i = 0; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		for (int j = 1; j < h - 2; j += 2)
		{
			for (int i = 1; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationOnOddPoints()
	{
#pragma omp parallel for
		for (int j = 1; j < h - 2; j += 2)
		{
			for (int i = 1; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationOnEvenPoints()
	{
#pragma omp parallel for
		for (int j = 0; j < h - 2; j += 2)
		{
			for (int i = 0; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationPostInterpolation()
	{
		GaussSeidelIterationOnSemiOddPoints();
		GaussSeidelIterationOnOddPoints();
		GaussSeidelIterationOnEvenPoints();
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIteration()
	{

		//Red Iteration
#pragma omp parallel for
		for (int j = 0; j < h - 2; j++)
		{
			for (int i = j % 2; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		//Black Iteration
		for (int j = 0; j < h - 2; j++)
		{
			for (int i = 1 - j % 2; i < w - 2; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

	virtual void SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity) = 0;

//-----------------------------------------------------------------------------

	void SolveWithGaussSeidel(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			F.Scale(-1);
			CalculateJ();
			
			U.SetAll(0);
			int gsItr = 0;
			while (true)
			{
				if (verbosity > 1)
					std::cout << "Gauss Seidel Itr: " << gsItr;
				ComputeResidualIntoO();
				double norm = O.Norm();
				if (verbosity > 1)
					std::cout << ", r: " << norm << std::endl;
				if (norm < maxNorm)
					break;
				GaussSeidelIteration();
				gsItr++;
				if (gsItr > maxItr)
					break;
			}

			if (gsItr > 0)
				std::cout << ", requiring " << gsItr << " GS iterations";
			double correctionNorm = U.Norm();
			if (verbosity > 0)
				std::cout << ", norm: " << correctionNorm << std::endl;
			if (correctionNorm < maxNorm)
				break;
			W.AddCorrectionToHaloedData(U, 1, minValue, maxValue);
			itr++;
		}
	}

//-----------------------------------------------------------------------------

	void SolveWithPARDISO(double maxNorm, double minValue, double maxValue, int verbosity)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			F.Scale(-1);
			CalculateJ();
			JS.LoadValues(J_cc, J_xp, J_xm, J_yp, J_ym);

			U.SetAll(0);
			SolveSparseSystemWithPARDISO(
				JS.matrixSize, 
				JS.GetValuesPtr(),
				JS.GetRowOffsetsPtr(),
				JS.GetColumnIndicesPtr(), 
				U.GetRawPtr(), 
				F.GetRawPtr(), 
				false);
			double correctionNorm = U.Norm();
			if (verbosity > 0)
				std::cout << ", norm: " << correctionNorm << std::endl;
			if (correctionNorm < maxNorm)
				break;
			W.AddCorrectionToHaloedData(U, 1, minValue, maxValue);
			itr++;
		}
	}

	//-----------------------------------------------------------------------------

	inline double GetW(int x, int y)
	{
		return W.Get(x, y);
	}
	inline double GetZ(int x, int y)
	{
		return GetB(x, y) + GetW(x, y);
	}
	inline double GetB(int x, int y)
	{
		return B.Get(x, y);
	}
	inline double GetR(int x, int y)
	{
		return R.Get(x, y);
	}
	inline void SetW(int x, int y, double value)
	{
		W.Set(x, y, value);
	}
	inline void SetB(int x, int y, double value)
	{
		B.Set(x, y, value);
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

//-----------------------------------------------------------------------------

	void RestrictRBWSimple(SolverDataBase<e - 1>& coarse)
	{
		for (int i = 0; i < coarse.w; i++)
		{
			for (int j = 0; j < coarse.h; j++)
			{
				coarse.SetB(i, j, GetB(2 * i, 2 * j));
				coarse.SetR(i, j, GetR(2 * i, 2 * j));
				coarse.SetW(i, j, GetW(2 * i, 2 * j));
			}
		}
	}

//-----------------------------------------------------------------------------

	void RestrictRBWTo(SolverDataBase<e - 1>& coarse)
	{
		for (int i = 0; i < coarse.w; i++)
		{
			for (int j = 0; j < coarse.h; j++)
			{
				if (i == 0 && j == 0)
				{
					coarse.SetB(i, j, GetB(2 * i, 2 * j));
					coarse.SetR(i, j, GetR(2 * i, 2 * j));
					coarse.SetW(i, j, GetW(2 * i, 2 * j));
				}
				else if (i == 0 && j == coarse.h - 1)
				{
					coarse.SetB(i, j, GetB(2 * i, 2 * j));
					coarse.SetR(i, j, GetR(2 * i, 2 * j));
					coarse.SetW(i, j, GetW(2 * i, 2 * j));
				}
				else if (i == coarse.w - 1 && j == 0)
				{
					coarse.SetB(i, j, GetB(2 * i, 2 * j));
					coarse.SetR(i, j, GetR(2 * i, 2 * j));
					coarse.SetW(i, j, GetW(2 * i, 2 * j));
				}
				else if (i == coarse.w - 1 && j == coarse.h - 1)
				{
					coarse.SetB(i, j, GetB(2 * i, 2 * j));
					coarse.SetR(i, j, GetR(2 * i, 2 * j));
					coarse.SetW(i, j, GetW(2 * i, 2 * j));
				}
				else if (i == 0)
				{
					double bVal = 0;
					bVal += GetB(2 * i, 2 * j - 1) / 4;
					bVal += GetB(2 * i, 2 * j) / 2;
					bVal += GetB(2 * i, 2 * j + 1) / 4;
					coarse.SetB(i, j, bVal);

					double rVal = 0;
					rVal += GetR(2 * i, 2 * j - 1) / 4;
					rVal += GetR(2 * i, 2 * j) / 2;
					rVal += GetR(2 * i, 2 * j + 1) / 4;
					coarse.SetR(i, j, rVal);

					double wVal = 0;
					wVal += GetB(2 * i, 2 * j - 1) / 4;
					wVal += GetB(2 * i, 2 * j) / 2;
					wVal += GetB(2 * i, 2 * j + 1) / 4;
					coarse.SetB(i, j, wVal);
				}
				else if (i == coarse.w - 1)
				{
					double bVal = 0;
					bVal += GetB(2 * i, 2 * j - 1) / 4;
					bVal += GetB(2 * i, 2 * j) / 2;
					bVal += GetB(2 * i, 2 * j + 1) / 4;
					coarse.SetB(i, j, bVal);

					double rVal = 0;
					rVal += GetR(2 * i, 2 * j - 1) / 4;
					rVal += GetR(2 * i, 2 * j) / 2;
					rVal += GetR(2 * i, 2 * j + 1) / 4;
					coarse.SetR(i, j, rVal);

					double wVal = 0;
					wVal += GetB(2 * i, 2 * j - 1) / 4;
					wVal += GetB(2 * i, 2 * j) / 2;
					wVal += GetB(2 * i, 2 * j + 1) / 4;
					coarse.SetB(i, j, wVal);
				}
				else if (j == 0)
				{
					double bVal = 0;
					bVal += GetB(2 * i - 1, 2 * j) / 4;
					bVal += GetB(2 * i, 2 * j) / 2;
					bVal += GetB(2 * i + 1, 2 * j) / 4;
					coarse.SetB(i, j, bVal);

					double rVal = 0;
					rVal += GetR(2 * i - 1, 2 * j) / 4;
					rVal += GetR(2 * i, 2 * j) / 2;
					rVal += GetR(2 * i + 1, 2 * j) / 4;
					coarse.SetR(i, j, rVal);

					double wVal = 0;
					wVal += GetW(2 * i - 1, 2 * j) / 4;
					wVal += GetW(2 * i, 2 * j) / 2;
					wVal += GetW(2 * i + 1, 2 * j) / 4;
					coarse.SetW(i, j, wVal);
				}
				else if (j == coarse.h - 1)
				{
					double bVal = 0;
					bVal += GetB(2 * i - 1, 2 * j) / 4;
					bVal += GetB(2 * i, 2 * j) / 2;
					bVal += GetB(2 * i + 1, 2 * j) / 4;
					coarse.SetB(i, j, bVal);

					double rVal = 0;
					rVal += GetR(2 * i - 1, 2 * j) / 4;
					rVal += GetR(2 * i, 2 * j) / 2;
					rVal += GetR(2 * i + 1, 2 * j) / 4;
					coarse.SetR(i, j, rVal);

					double wVal = 0;
					wVal += GetW(2 * i - 1, 2 * j) / 4;
					wVal += GetW(2 * i, 2 * j) / 2;
					wVal += GetW(2 * i + 1, 2 * j) / 4;
					coarse.SetW(i, j, wVal);
				}
				else
				{
					double bVal = 0;
					bVal += GetB(2 * i - 1, 2 * j - 1) / 16;
					bVal += GetB(2 * i, 2 * j - 1) / 8;
					bVal += GetB(2 * i + 1, 2 * j - 1) / 16;
					bVal += GetB(2 * i - 1, 2 * j) / 8;
					bVal += GetB(2 * i, 2 * j) / 4;
					bVal += GetB(2 * i + 1, 2 * j) / 8;
					bVal += GetB(2 * i - 1, 2 * j + 1) / 16;
					bVal += GetB(2 * i, 2 * j + 1) / 8;
					bVal += GetB(2 * i + 1, 2 * j + 1) / 16;
					coarse.SetB(i, j, bVal);

					double rVal = 0;
					rVal += GetR(2 * i - 1, 2 * j - 1) / 16;
					rVal += GetR(2 * i, 2 * j - 1) / 8;
					rVal += GetR(2 * i + 1, 2 * j - 1) / 16;
					rVal += GetR(2 * i - 1, 2 * j) / 8;
					rVal += GetR(2 * i, 2 * j) / 4;
					rVal += GetR(2 * i + 1, 2 * j) / 8;
					rVal += GetR(2 * i - 1, 2 * j + 1) / 16;
					rVal += GetR(2 * i, 2 * j + 1) / 8;
					rVal += GetR(2 * i + 1, 2 * j + 1) / 16;
					coarse.SetR(i, j, rVal);

					double wVal = 0;
					wVal += GetW(2 * i - 1, 2 * j - 1) / 16;
					wVal += GetW(2 * i, 2 * j - 1) / 8;
					wVal += GetW(2 * i + 1, 2 * j - 1) / 16;
					wVal += GetW(2 * i - 1, 2 * j) / 8;
					wVal += GetW(2 * i, 2 * j) / 4;
					wVal += GetW(2 * i + 1, 2 * j) / 8;
					wVal += GetW(2 * i - 1, 2 * j + 1) / 16;
					wVal += GetW(2 * i, 2 * j + 1) / 8;
					wVal += GetW(2 * i + 1, 2 * j + 1) / 16;
					coarse.SetW(i, j, wVal);
				}
			}
		}
	}

//-----------------------------------------------------------------------------

	void InterpolateWFrom(SolverDataBase<e - 1>& coarse)
	{
		for (int j = 1; j < h - 1; j++)
		{
			for (int i = 1; i < w - 1; i++)
			{
				double value = 0;
				int ci = i / 2;
				int cj = j / 2;
				if( i % 2 == 0 && j % 2 == 0)
				{
					value = coarse.GetW(ci, cj);
					SetW(i, j, value);
				}
				else if (i % 2 == 0)
				{
					value += coarse.GetW(ci, cj) * 0.5;
					value += coarse.GetW(ci, cj + 1) * 0.5;
					SetW(i, j, value);
				}
				else if (j % 2 == 0)
				{
					value += coarse.GetW(ci, cj) * 0.5;
					value += coarse.GetW(ci + 1, cj) * 0.5;
					SetW(i, j, value);
				}
				else
				{
					value += coarse.GetW(ci, cj) * 0.25;
					value += coarse.GetW(ci + 1, cj) * 0.25;
					value += coarse.GetW(ci, cj + 1) * 0.25;
					value += coarse.GetW(ci + 1, cj + 1) * 0.25;
					SetW(i, j, value);
				}
			}
		}
	}

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
	DataPointer<w, h> B; //Bed Height and boundary condition vectors
	DataPointer<w, h> W; //Water Height and boundary condition vectors
	DataPointer<w, h> R; //Rainfall Values

	DataPointer<w - 2, h - 2> F; //buffer for evaluations of F_ij(W)
	DataPointer<w - 2, h - 2> U; //buffer for correction that we calculate
	DataPointer<w - 2, h - 2> O; //buffer for random stuff - residuals, for instance

	/*
	* The Jacobian of F is a pentadiagonal matrix with each diagonal
	* corresponding to one spoke of the 5 point stencil.
	* We store each diagonal (each spoke of the stencil) in a separate
	* DataPointer for maximum efficiency in sparsity
	*
	*/
	DataPointer<w - 2, h - 2> J_cc; //Center
	DataPointer<w - 2, h - 2> J_xp; //i+1, j
	DataPointer<w - 2, h - 2> J_xm; //i-1, j
	DataPointer<w - 2, h - 2> J_yp; //i, j+1
	DataPointer<w - 2, h - 2> J_ym; //i, j-1

	FiveStencilCSRMatrix<w - 2, h - 2> JS;
};

template <int e>
class SolverData : public SolverDataBase<e>
{
public:
	SolverData(double gravity, double drag, double spacing) :
		SolverDataBase<e>(gravity, drag, spacing),
		child(gravity, drag, spacing * 2)
	{

	}
	void SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		this->RestrictRBWSimple(child);
		child.SolveWithPseudoMultigrid(maxNorm, minValue, maxValue, maxItr, verbosity);
		this->InterpolateWFrom(child);
		this->GaussSeidelIterationPostInterpolation();

		std::cout << "Solving on Mesh level " << e << std::endl;
		this->SolveWithGaussSeidel(maxNorm, minValue, maxValue, maxItr, verbosity);
	}
private:
	SolverData<e - 1> child;
};

template <>
class SolverData<2> : public SolverDataBase<2>
{
public:
	SolverData(double gravity, double drag, double spacing) :
		SolverDataBase<2>(gravity, drag, spacing)	
	{
	}
	void SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		std::cout << "Solving on Mesh level " << 2 << std::endl;
		this->SolveWithGaussSeidel(maxNorm, minValue, maxValue, maxItr, verbosity);
	}
};