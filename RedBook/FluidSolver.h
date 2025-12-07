#pragma once

#include "Heightmap.h"
#include "Depthmap.h"
#include <filesystem>

template <int e>
class SolverData;

class FluidSolver
{
public:
	FluidSolver(int LOD);
	~FluidSolver();
	void LoadDefaultRainfallData();
	void LoadDefaultDepthGuess();
	void LoadHeightmap(RegionalDataLoc where, HeightmapManager& source);
	void LoadDepthmap(RegionalDataLoc where, DepthmapManager& source, bool demand);
	void VisualizeData();
	void Solve();
	void FullSolutionCycle(RegionalDataLoc where, HeightmapManager& heights, DepthmapManager& depths);
	void SaveDepthmap(RegionalDataLoc where, DepthmapManager& source);
	void SmoothHeightmap(float smoothingFactor);
	void RecursiveSolutionCycle(RegionalDataLoc where, HeightmapManager& heights, DepthmapManager& depths, int recursiveDepth);

private:
	void Visualize(std::filesystem::path outPath, std::string name, int mode);
public:
	static constexpr int dim = 8;
	static constexpr double GRAVITY = 9.8;
	static constexpr double DRAG = 0.7;
	static constexpr double ANNUAL_RAINFALL = 0.000000049; //meters per second, about 64 inches / year
private:
	SolverData<dim>* data;
};

