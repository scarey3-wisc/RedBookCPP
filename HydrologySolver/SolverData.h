#pragma once
#include "OperationModes.h"


#include "MatrixOps.h"
#include "MultigridOps.h"
#include "CSRMatrix.h"

#ifdef USE_PARDISO
#include "mkl.h"
#endif

#include <iostream>
#include <omp.h>
#include <limits>
#include "FileGlobals.h"

template <int e>
class SolverDataBase
{
public:
	SolverDataBase(double gravity, double drag, double sideLength) :
		g(gravity), l(drag), s(sideLength / (wI - 1)), fc(gravity / drag), s2(s*s),
		B(DataPointer<wO, hO>::AllocateNewPointer()),
		W(DataPointer<wO, hO>::AllocateNewPointer()),
		R(DataPointer<wO, hO>::AllocateNewPointer()),
		F(DataPointer<wI, hI>::AllocateNewPointer()),
		U(DataPointer<wI, hI>::AllocateNewPointer()),
		O(DataPointer<wI, hI>::AllocateNewPointer()),
		J_cc(DataPointer<wI, hI>::AllocateNewPointer()),
		J_xp(DataPointer<wI, hI>::AllocateNewPointer()),
		J_xm(DataPointer<wI, hI>::AllocateNewPointer()),
		J_yp(DataPointer<wI, hI>::AllocateNewPointer()),
		J_ym(DataPointer<wI, hI>::AllocateNewPointer())
#ifdef USE_PARDISO
		,
		JS(FiveStencilCSRMatrix<wI, hI>::AllocateNewPointer()),
		pardisoInit(false)
#endif
	{

	}

//-----------------------------------------------------------------------------

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
#ifdef USE_PARDISO
		JS.Deallocate();
#endif
	}

	// Disable copying:
	SolverDataBase(const SolverDataBase&) = delete;
	SolverDataBase& operator=(const SolverDataBase&) = delete;

//-----------------------------------------------------------------------------

	//outer and inner dimensions of the underlying 2D grid.
	//there's a single layer of data nodes whose values
	//we take as fixed boundary conditions.
	static constexpr int wO = (1 << e) + 1;
	static constexpr int hO = (1 << e) + 1;
	static constexpr int wI = (1 << e) - 1;
	static constexpr int hI = (1 << e) - 1;

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
	inline void SetSpacing(double domainSideLength)
	{
		s = domainSideLength / (wI - 1);
		s2 = s * s;
	}
	inline void SetGravity(double gravity)
	{
		g = gravity;
		fc = g / l;
	}
	inline void SetDrag(double drag)
	{
		l = drag;
		fc = g / l;
	}

//-----------------------------------------------------------------------------

	void ComputeResidualIntoO()
	{
#pragma omp parallel for
		for (int j = 0; j < hI; j++)
		{
			for (int i = 0; i < wI; i++)
			{
				double result = J_cc.Get(i, j) * U.Get(i, j);
				if (i > 0)
					result += J_xm.Get(i, j) * U.Get(i - 1, j);
				if (i < wI - 1)
					result += J_xp.Get(i, j) * U.Get(i + 1, j);
				if (j > 0)
					result += J_ym.Get(i, j) * U.Get(i, j - 1);
				if (j < hI - 1)
					result += J_yp.Get(i, j) * U.Get(i, j + 1);

				result -= F.Get(i, j);
				result *= -1;
				O.Set(i, j, result);
			}
		}
	}

//-----------------------------------------------------------------------------

	void CalculateF()
	{
#pragma omp parallel for
		for (int j = 1; j < hO - 1; j++)
		{
			for (int i = 1; i < wO - 1; i++)
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
		for (int j = 1; j < hO - 1; j++)
		{
			for (int i = 1; i < wO - 1; i++)
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

#ifdef USE_PRESMOOTHED_PARDISO
	virtual int SolveWithPreSmoothedPARDISO(double maxNorm, double minValue, double maxValue, int verbosity) = 0;
#endif

#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
	virtual int SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity) = 0;
#endif

#ifdef USE_PRESMOOTHED_MULTIGRID
	virtual int SolveWithMultigrid(double maxNorm, double minValue, double maxValue, int maxVCycles, int maxGSItr, int verbosity) = 0;

	virtual void PrepareMultigridJacobians() = 0;

	virtual void VCycle(double maxGSNorm, int maxGSItr, int gsItrDown, int gsItrUp) = 0;
#endif

//-----------------------------------------------------------------------------

	int SolveWithGaussSeidel(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			CalculateJ();
			F.Scale(-1);

			double residualNorm = F.Norm();
			if (verbosity > 0)
				std::cout << ", residual norm: " << residualNorm;
			if (residualNorm < maxNorm)
			{
				if(verbosity > 0)
					std::cout << std::endl;
				break;
			}
			
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

			if (gsItr > 0 && verbosity > 0)
				std::cout << ", requiring " << gsItr << " GS iterations" << std::endl;
			W.AddCorrectionToHaloedData(U, 1, minValue, maxValue);
			itr++;
			if (itr > 2000)
			{
				std::cout << "Failed to converge after 2000 Newton Iterations" << std::endl;
				break;
			}
		}
		return itr;
	}

//-----------------------------------------------------------------------------

#ifdef USE_PARDISO
	int SolveWithPARDISO(double maxNorm, double minValue, double maxValue, int verbosity)
	{
		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			CalculateF();
			CalculateJ();
			J_cc.AddAbs(F, -2.0);
			JS.LoadValues(J_cc, J_xp, J_xm, J_yp, J_ym);
			F.Scale(-1);

			double residualNorm = F.Norm();
			if (verbosity > 0)
				std::cout << ", residual norm: " << residualNorm << std::endl;
			if (residualNorm < maxNorm)
			{
				if (verbosity > 0)
					std::cout << std::endl;
				break;
			}

			U.SetAll(0);
			SolveSparseSystemWithPARDISO(
				JS.matrixSize, 
				JS.GetValuesPtr(),
				JS.GetRowOffsetsPtr(),
				JS.GetColumnIndicesPtr(), 
				U.GetRawPtr(), 
				F.GetRawPtr(), 
				false);
			W.AddCorrectionToHaloedData(U, 1, minValue, maxValue);
			itr++;
			if (itr > 2000)
			{
				std::cout << "Failed to converge after 2000 Newton Iterations" << std::endl;
				break;
			}
		}
		TerminatePARDISOStructures(JS.matrixSize, JS.GetRowOffsetsPtr(), JS.GetColumnIndicesPtr(), false);
		return itr;
	}
#endif	

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
		for (int j = 0; j < hI; j++)
		{
			for (int i = j % 2; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		//Black Iteration
		for (int j = 0; j < hI; j++)
		{
			for (int i = 1 - j % 2; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

private:

#ifdef USE_PARDISO

//-----------------------------------------------------------------------------

	void OnFirstPARDISOSolveAttempt(int matrixSize, double* vals, int* rowOffsets, int* colIndices, bool verbose)
	{

		MKL_INT n = matrixSize; // Matrix size
		MKL_INT mtype = 1;        // Real and structurally symmetric matrix
		MKL_INT nrhs = 1;         // Number of right hand sides
		// should be "int" when using 32-bit architectures, or "long int"
		// for 64-bit architectures. void* should be OK in both cases
		MKL_INT maxfct, mnum, phase, error, msglvl;
		MKL_INT i;
		float ddum;               // Scalar dummy (PARDISO needs it)
		MKL_INT idum;             // Integer dummy (PARDISO needs it)

		maxfct = 1;           // Maximum number of numerical factorizations.
		mnum = 1;             // Which factorization to use.
		msglvl = verbose ? 1 : 0;           // Print statistical information in file
		error = 0;            // Initialize error flag

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
		iparm[20] = 1;		  // USe some pivoting to help with singularity
		iparm[23] = 1;        // Two-level factorization*/
		iparm[26] = 1;        // Check matrix for errors
		iparm[27] = 0;        // Use double precision
		iparm[34] = 1;        // Use zero-based indexing

		// Initialize the internal solver memory pointer. This is only
		// necessary for the FIRST call of the PARDISO solver
		for (i = 0; i < 64; i++)
		{
			pardisoPT[i] = 0;
		}

		// Reordering and Symbolic Factorization. This step also allocates
		// all memory that is necessary for the factorization
		phase = 11;
		PARDISO(pardisoPT, &maxfct, &mnum, &mtype, &phase, &n,
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

		pardisoInit = true;
	}

//-----------------------------------------------------------------------------

	void SolveSparseSystemWithPARDISO(int matrixSize, double* vals, int* rowOffsets, int* colIndices, double* x, double* f, bool verbose)
	{
		MKL_INT n = matrixSize; // Matrix size
		MKL_INT mtype = 1;        // Real and structurally symmetric matrix
		MKL_INT nrhs = 1;         // Number of right hand sides
		// should be "int" when using 32-bit architectures, or "long int"
		// for 64-bit architectures. void* should be OK in both cases
		MKL_INT maxfct, mnum, phase, error, msglvl;
		float ddum;               // Scalar dummy (PARDISO needs it)
		MKL_INT idum;             // Integer dummy (PARDISO needs it)

		maxfct = 1;           // Maximum number of numerical factorizations.
		mnum = 1;             // Which factorization to use.
		msglvl = verbose ? 1 : 0;           // Print statistical information in file
		error = 0;            // Initialize error flag

		if (!pardisoInit)
			OnFirstPARDISOSolveAttempt(matrixSize, vals, rowOffsets, colIndices, verbose);

		// Numerical factorization
		phase = 22;
		PARDISO(pardisoPT, &maxfct, &mnum, &mtype, &phase, &n,
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
		PARDISO(pardisoPT, &maxfct, &mnum, &mtype, &phase, &n,
			vals, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, static_cast<void*>(f), x, &error);
		if (error != 0)
			throw std::runtime_error("PARDISO error during solution phase");

		if (verbose)
		{
			std::cout << "Solve completed ... " << std::endl;
		}
	}

//-----------------------------------------------------------------------------

	void TerminatePARDISOStructures(int matrixSize, int* rowOffsets, int* colIndices, bool verbose)
	{

		MKL_INT n = matrixSize; // Matrix size
		MKL_INT mtype = 1;        // Real and structurally symmetric matrix
		MKL_INT nrhs = 1;         // Number of right hand sides
		// should be "int" when using 32-bit architectures, or "long int"
		// for 64-bit architectures. void* should be OK in both cases
		MKL_INT maxfct, mnum, phase, error, msglvl;
		float ddum;               // Scalar dummy (PARDISO needs it)
		MKL_INT idum;             // Integer dummy (PARDISO needs it)

		maxfct = 1;           // Maximum number of numerical factorizations.
		mnum = 1;             // Which factorization to use.
		msglvl = verbose ? 1 : 0;           // Print statistical information in file
		error = 0;            // Initialize error flag
		// Termination and release of memory.
		phase = -1;           // Release internal memory
		PARDISO(pardisoPT, &maxfct, &mnum, &mtype, &phase, &n,
			&ddum, rowOffsets, colIndices,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

		pardisoInit = false;
	}

#endif

//-----------------------------------------------------------------------------

	void SmoothIndex(int i, int j)
	{
		double rhs = F.Get(i, j);
		if (i > 0)
			rhs -= J_xm.Get(i, j) * U.Get(i - 1, j);
		if (i < wI - 1)
			rhs -= J_xp.Get(i, j) * U.Get(i + 1, j);
		if (j > 0)
			rhs -= J_ym.Get(i, j) * U.Get(i, j - 1);
		if (j < hI - 1)
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
		for (int j = 0; j < hI; j += 2)
		{
			for (int i = 1; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		for (int j = 1; j < hI; j += 2)
		{
			for (int i = 0; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
#pragma omp parallel for
		for (int j = 1; j < hI; j += 2)
		{
			for (int i = 1; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationOnOddPoints()
	{
#pragma omp parallel for
		for (int j = 1; j < hI; j += 2)
		{
			for (int i = 1; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	void GaussSeidelIterationOnEvenPoints()
	{
#pragma omp parallel for
		for (int j = 0; j < hI; j += 2)
		{
			for (int i = 0; i < wI; i += 2)
			{
				SmoothIndex(i, j);
			}
		}
	}

//-----------------------------------------------------------------------------

	double g; // gravitational constant
	double l; // drag constant
	double s; // spacing constant

	double fc; // the quantity g/2l, an important flux constant
	double s2; // spacing squared is often an important quantity

#ifdef USE_PARDISO
	void* pardisoPT[64];	// Internal solver memory pointer pt for PARDISO
	MKL_INT iparm[64];        // Pardiso control parameters
	bool pardisoInit;
#endif
	
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

public:
	DataPointer<wO, hO> B; //Bed Height and boundary condition vectors
	DataPointer<wO, hO> W; //Water Height and boundary condition vectors
	DataPointer<wO, hO> R; //Rainfall Values

	DataPointer<wI, hI> F; //buffer for evaluations of F_ij(W)
	DataPointer<wI, hI> U; //buffer for correction that we calculate
	DataPointer<wI, hI> O; //buffer for random stuff - residuals, for instance

	/*
	* The Jacobian of F is a pentadiagonal matrix with each diagonal
	* corresponding to one spoke of the 5 point stencil.
	* We store each diagonal (each spoke of the stencil) in a separate
	* DataPointer for maximum efficiency in sparsity
	*/
	DataPointer<wI, hI> J_cc; //Center
	DataPointer<wI, hI> J_xp; //i+1, j
	DataPointer<wI, hI> J_xm; //i-1, j
	DataPointer<wI, hI> J_yp; //i, j+1
	DataPointer<wI, hI> J_ym; //i, j-1

#ifdef USE_PARDISO
	FiveStencilCSRMatrix<wI, hI> JS;
#endif
};

template <int e>
class SolverData : public SolverDataBase<e>
{
public:
	SolverData(double gravity, double drag, double sideLength) :
		SolverDataBase<e>(gravity, drag, sideLength)
#ifdef RECURSIVE_STRUCTURE_REQUIRED
		,
		child(gravity, drag, sideLength)
#endif
	{

	}
#ifdef RECURSIVE_STRUCTURE_REQUIRED
	inline void SetSpacing(double domainSideLength)
	{
		SolverDataBase<e>::SetSpacing(domainSideLength);
		child.SetSpacing(domainSideLength);
	}
	inline void SetGravity(double gravity)
	{
		SolverDataBase<e>::SetGravity(gravity);
		child.SetGravity(gravity);
	}
	inline void SetDrag(double drag)
	{
		SolverDataBase<e>::SetDrag(drag);
		child.SetDrag(drag);
	}
#endif

#ifdef USE_PRESMOOTHED_PARDISO
	int SolveWithPreSmoothedPARDISO(double maxNorm, double minValue, double maxValue, int verbosity)
	{
		MultigridOps::Restrict<e>(this->R, child.R);
		MultigridOps::Restrict<e>(this->B, child.B);
		MultigridOps::Restrict<e>(this->W, child.W);
		child.SolveWithPreSmoothedPARDISO(maxNorm, minValue, maxValue, verbosity);
		MultigridOps::Interpolate<e>(this->W, child.W);
		this->CalculateJ();
		this->CalculateF();
		this->GaussSeidelIterationPostInterpolation();

		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << e << std::endl;
		return this->SolveWithPARDISO(maxNorm, minValue, maxValue, verbosity);
	}
#endif

#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
	int SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		MultigridOps::Restrict<e>(this->R, child.R);
		MultigridOps::Restrict<e>(this->B, child.B);
		MultigridOps::Restrict<e>(this->W, child.W);
		child.SolveWithPseudoMultigrid(maxNorm, minValue, maxValue, maxItr, verbosity);
		MultigridOps::Interpolate<e>(this->W, child.W);
		this->CalculateJ();
		this->CalculateF();
		this->GaussSeidelIterationPostInterpolation();

		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << e << std::endl;
		return this->SolveWithGaussSeidel(maxNorm, minValue, maxValue, maxItr, verbosity);
	}
#endif

//-----------------------------------------------------------------------------

#ifdef USE_PRESMOOTHED_MULTIGRID
	int SolveWithMultigrid(double maxNorm, double minValue, double maxValue, int maxVCycles, int maxGSItr, int verbosity)
	{
		MultigridOps::Restrict<e>(this->R, child.R);
		MultigridOps::Restrict<e>(this->B, child.B);
		MultigridOps::Restrict<e>(this->W, child.W);
		child.SolveWithMultigrid(maxNorm, minValue, maxValue, maxVCycles, maxGSItr, verbosity);
		MultigridOps::Interpolate<e>(this->W, child.W);
		this->CalculateJ();
		this->CalculateF();
		//this->GaussSeidelIterationPostInterpolation();

		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << e << std::endl;

		int itr = 0;
		while (true)
		{
			if (verbosity > 0)
				std::cout << "Newton Itr " << itr;
			
			PrepareMultigridJacobians();
			this->CalculateF();
			this->F.Scale(-1);

			double residualNorm = this->F.Norm();
			if (verbosity > 0)
				std::cout << ", residual norm: " << residualNorm;
			if (residualNorm < maxNorm && itr > 0)
			{
				if(verbosity > 0)
					std::cout << std::endl;
				break;
			}

			this->U.SetAll(0);

			int vItr = 0;
			while (true)
			{
				if (verbosity > 1)
					std::cout << "VCycle: " << vItr;
				this->ComputeResidualIntoO();
				double norm = this->O.Norm();
				if (verbosity > 1)
					std::cout << ", r: " << norm << std::endl;
				if (norm < maxNorm && vItr > 0)
					break;
				VCycle(maxNorm, maxGSItr, 4, 4);
				vItr++;
				if (vItr >= maxVCycles)
					break;
			}
			if (vItr > 0 && verbosity > 0)
				std::cout << ", requiring " << vItr << " V Cycles" << std::endl;
			this->W.AddCorrectionToHaloedData(this->U, 1, minValue, maxValue);
			itr++;
			if (itr > 2000)
			{
				std::cout << "Failed to converge after 2000 Newton Iterations" << std::endl;
				break;
			}
		}
		return itr;
	}
	void PrepareMultigridJacobians()
	{
		this->CalculateJ();
		MultigridOps::Restrict<e>(this->R, child.R);
		MultigridOps::Restrict<e>(this->B, child.B);
		MultigridOps::Restrict<e>(this->W, child.W);
		child.PrepareMultigridJacobians();
	}
	void VCycle(double maxGSNorm, int maxGSItr, int gsItrDown, int gsItrUp)
	{
		this->ComputeResidualIntoO();
		double currNorm = this->O.Norm();
		for (int i = 0; i < gsItrDown; i++)
		{
			this->GaussSeidelIteration();
		}
		this->ComputeResidualIntoO();
		currNorm = this->O.Norm();
		MultigridOps::Restrict<e>(this->O, child.F);
		child.U.SetAll(0);
		child.VCycle(maxGSNorm, maxGSItr, gsItrDown, gsItrUp);
		MultigridOps::Interpolate<e>(this->O, child.U);
		this->U.Add(this->O, 1.0);

		this->ComputeResidualIntoO();
		currNorm = this->O.Norm();
		for(int i = 0; i < gsItrUp; i++)
		{
			this->GaussSeidelIteration();
		}

		this->ComputeResidualIntoO();
		currNorm = this->O.Norm();
	}
#endif

#ifdef RECURSIVE_STRUCTURE_REQUIRED
private:
	SolverData<e - 1> child;
#endif

};

#ifdef RECURSIVE_STRUCTURE_REQUIRED

template <>
class SolverData<4> : public SolverDataBase<4>
{
public:
	SolverData(double gravity, double drag, double sideLength) :
		SolverDataBase<4>(gravity, drag, sideLength)	
	{

	}

#ifdef USE_PRESMOOTHED_PARDISO
	int SolveWithPreSmoothedPARDISO(double maxNorm, double minValue, double maxValue, int verbosity)
	{
		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << 4 << std::endl;
		return this->SolveWithPARDISO(maxNorm, minValue, maxValue, verbosity);
	}
#endif

#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
	int SolveWithPseudoMultigrid(double maxNorm, double minValue, double maxValue, int maxItr, int verbosity)
	{
		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << 4 << std::endl;
		return this->SolveWithGaussSeidel(maxNorm, minValue, maxValue, maxItr, verbosity);
	}
#endif

#ifdef USE_PRESMOOTHED_MULTIGRID
	int SolveWithMultigrid(double maxNorm, double minValue, double maxValue, int maxVCycles, int maxGSItr, int verbosity)
	{
		if(verbosity > 0)
			std::cout << "Solving on Mesh level " << 4 << std::endl;
		return this->SolveWithGaussSeidel(maxNorm, minValue, maxValue, maxGSItr, verbosity);
	}
	void PrepareMultigridJacobians()
	{
		this->CalculateJ();
	}
	void VCycle(double maxGSNorm, int maxGSItr, int gsItrDown, int gsItrUp)
	{
		int gsItr = 0;
		for (int i = 0; i < maxGSItr; i++)
		{
			ComputeResidualIntoO();
			double norm = O.Norm();
			if (norm < maxGSNorm)
				break;
			GaussSeidelIteration();
			gsItr++;
			if (gsItr > maxGSItr)
				break;
		}
	}
#endif
};
#endif