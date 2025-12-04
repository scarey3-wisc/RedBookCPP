#pragma once
#include "../HydrologySolver/SolverData.h"
#include "Heightmap.h"

class FluidSolver
{
public:
	void LoadDefaultRainfallData();
	void LoadDefaultDepthGuess();
public:
	static constexpr int dim = 8;
	static constexpr double GRAVITY = 9.8;
	static constexpr double DRAG = 0.7;
	static constexpr double ANNUAL_RAINFALL = 0.000000049; //meters per second, about 64 inches / year
private:
	SolverData<dim> data;
};

