#include "FluidSolver.h"
#include "../HydrologySolver/SolverData.h"
#include "RegionalMap.h"

FluidSolver::FluidSolver(int LOD)
{
    double meterDim = 1.0 * RegionalMap::METER_DIM;
    meterDim /= (1 << (LOD - 1));
    data = new SolverData<dim>(GRAVITY, DRAG, meterDim);
}
FluidSolver::~FluidSolver()
{
    if (data != nullptr)
    {
        delete data;
        data = nullptr;
    }
}
void 
FluidSolver::FullSolutionCycle(RegionalDataLoc where, HeightmapManager& heights, DepthmapManager& depths)
{
    double meterDim = 1.0 * RegionalMap::METER_DIM;
    meterDim /= (1.0 * (1 << (where.LOD - 1)));
    data->SetSpacing(meterDim);
    LoadDefaultDepthGuess();
    LoadDefaultRainfallData();
    LoadHeightmap(where, heights);
    LoadDepthmap(where, depths, false);
    VisualizeData();
    Solve();
    VisualizeData();
    SaveDepthmap(where, depths);
}

void
FluidSolver::RecursiveSolutionCycle(RegionalDataLoc where, HeightmapManager& heights, DepthmapManager& depths, int recursiveDepth)
{
    int wX = (int)where.x;
    int wY = (int)where.y;
    int wD = (int)where.LOD;
    std::cout << std::endl;
    std::cout << "||**************************************||" << std::endl;
    std::cout << "||****** SOLVING " << wX << ", " << wY << ", AT DEPTH " << wD << " ******|| " << std::endl;
    std::cout << "||**************************************||" << std::endl << std::endl;

    double meterDim = 1.0 * RegionalMap::METER_DIM;
    meterDim /= (1.0 * (1 << (where.LOD - 1)));
    data->SetSpacing(meterDim);
    LoadDefaultDepthGuess();
    LoadDefaultRainfallData();
    LoadHeightmap(where, heights);
    LoadDepthmap(where, depths, true);
    Solve();
    SaveDepthmap(where, depths);

    if (recursiveDepth <= 1)
        return;
    for (RegionalDataLoc child : where.GetHigherLODs())
    {
        RecursiveSolutionCycle(child, heights, depths, recursiveDepth - 1);
    }
}


void 
FluidSolver::LoadDefaultRainfallData()
{
    for (int j = 0; j < data->hO; j++)
    {
        for (int i = 0; i < data->wO; i++)
        {
            data->SetR(i, j, 0);
        }
    }
    for (int j = 1; j < data->hO - 1; j++)
    {
        for (int i = 1; i < data->wO - 1; i++)
        {
            data->SetR(i, j, 200 * ANNUAL_RAINFALL);
        }
    }
}
void 
FluidSolver::LoadDefaultDepthGuess()
{
    for (int j = 0; j < data->hO; j++)
    {
        for (int i = 0; i < data->wO; i++)
        {
            data->SetW(i, j, 0.0);
        }
    }
    for (int j = 1; j < data->hO - 1; j++)
    {
        for (int i = 1; i < data->wO - 1; i++)
        {
            data->SetW(i, j, 1.0);
        }
    }
}
void 
FluidSolver::SmoothHeightmap(float smoothingFactor)
{
    for (int j = 1; j < data->hO - 1; j++)
    {
        for (int i = 1 + j % 2; i < data->wO - 1; i+=2)
        {
            double laplacian = 0;
            laplacian += 4 * data->GetB(i, j);
            laplacian -= data->GetB(i - 1, j);
            laplacian -= data->GetB(i + 1, j);
            laplacian -= data->GetB(i, j - 1);
            laplacian -= data->GetB(i, j + 1);
            laplacian *= smoothingFactor;

            double newHeight = data->GetB(i, j) - laplacian;
            data->SetB(i, j, newHeight);
        }
    }

    for (int j = 1; j < data->hO - 1; j++)
    {
        for (int i = 2 - j % 2; i < data->wO - 1; i += 2)
        {
            double laplacian = 0;
            laplacian += 4 * data->GetB(i, j);
            laplacian -= data->GetB(i - 1, j);
            laplacian -= data->GetB(i + 1, j);
            laplacian -= data->GetB(i, j - 1);
            laplacian -= data->GetB(i, j + 1);
            laplacian *= smoothingFactor;

            double newHeight = data->GetB(i, j) - laplacian;
            data->SetB(i, j, newHeight);
        }
    }
}

void 
FluidSolver::SaveDepthmap(RegionalDataLoc where, DepthmapManager& source)
{
    Depthmap target = source.DemandRaster(where);
    for (int j = 0; j < data->hO; j++)
    {
        for (int i = 0; i < data->wO; i++)
        {
            double depth = data->GetW(i, j);
            target.SetDepth(i + 1, j + 1, (float)depth);
        }
    }
}

void 
FluidSolver::LoadDepthmap(RegionalDataLoc where, DepthmapManager& source, bool demand)
{
    double meterDim = 1.0 * RegionalMap::METER_DIM;
    meterDim /= where.GetNumSectionsPerSide();
    data->SetSpacing(meterDim);

    if (!source.DataAvailable(where, false) && !demand)
        return;

    Depthmap target = source.DemandRaster(where);
    for (int j = 0; j < data->hO; j++)
    {
        for (int i = 0; i < data->wO; i++)
        {
            float depth = target.GetDepth(i + 1, j + 1);
            data->SetW(i, j, (double) depth);
        }
    }
}

void 
FluidSolver::LoadHeightmap(RegionalDataLoc where, HeightmapManager& source)
{
    double meterDim = 1.0 * RegionalMap::METER_DIM;
    meterDim /= where.GetNumSectionsPerSide();
    data->SetSpacing(meterDim);
    Heightmap target = source.DemandRaster(where);
    for (int j = 0; j < data->hO; j++)
    {
        for (int i = 0; i < data->wO; i++)
        {
            double elev = target.GetElevation(i + 1, j + 1);
            if (elev == 0)
                elev = -100;
            data->SetB(i, j, elev);
        }
    }
}
void 
FluidSolver::VisualizeData()
{
    std::filesystem::path outPath = MY_PATH;
    Visualize(outPath, "waterways", 0);
    Visualize(outPath, "terrain", 1);
    Visualize(outPath, "surface", 2);
}
void 
FluidSolver::Visualize(std::filesystem::path outPath, std::string name, int mode)
{
    outPath.replace_filename(name + ".ppm");
    std::ofstream out(outPath.c_str());
    out << "P3\n" << data->wO << " " << data->hO << "\n255\n";
    double minB = std::numeric_limits<double>::max();
    double maxB = std::numeric_limits<double>::min();
    double minW = std::numeric_limits<double>::max();
    double maxW = std::numeric_limits<double>::min();
    for (int y = 0; y < data->hO; ++y) {
        for (int x = 0; x < data->wO; ++x) {
            if (data->GetW(x, y) < minW)
                minW = data->GetW(x, y);
            if (data->GetW(x, y) > maxW)
                maxW = data->GetW(x, y);
            if (mode == 2)
            {
                if (data->GetZ(x, y) < minB)
                    minB = data->GetZ(x, y);
                if (data->GetZ(x, y) > maxB)
                    maxB = data->GetZ(x, y);
            }
            else
            {
                if (data->GetB(x, y) < minB)
                    minB = data->GetB(x, y);
                if (data->GetB(x, y) > maxB)
                    maxB = data->GetB(x, y);
            }

        }
    }
    minW = 0;
    std::cout << "B range: " << minB << " - " << maxB << ", W range: " << minW << " - " << maxW << std::endl;
    for (int y = 0; y < data->hO; ++y) {
        for (int x = 0; x < data->wO; ++x) {
            double w = (data->GetW(x, y) - minW) / (maxW - minW);
            double t = 0;
            if (maxB > minB)
                t = (data->GetB(x, y) - minB) / (maxB - minB);
            if (mode == 2 && maxB > minB)
                t = (data->GetZ(x, y) - minB) / (maxB - minB);
            if (mode == 0)
                w = std::pow(w, 0.5);
            else
                w = 0;

            if (data->GetW(x, y) >= 5.0)
                w = 1;

            int wR = 0;
            int wG = 0;
            int wB = 255;
            int lR = 0;
            int lG = 255;
            int lB = 0;
            int hR = 255;
            int hG = 0;
            int hB = 0;
            if (mode > 0)
            {
                lR = 0;
                lG = 0;
                lB = 0;
                hR = 255;
                hG = 255;
                hB = 255;
            }

            int tR = (int)(lR * (1 - t) + hR * t);
            int tG = (int)(lG * (1 - t) + hG * t);
            int tB = (int)(lB * (1 - t) + hB * t);

            int r = (int)(tR * (1 - w) + wR * w);
            int g = (int)(tG * (1 - w) + wG * w);
            int b = (int)(tB * (1 - w) + wB * w);
            out << r << " " << g << " " << b << " ";
        }
        out << "\n";
    }
}
void
FluidSolver::Solve()
{
    if (data == nullptr)
        return;

#ifdef USE_PRESMOOTHED_PARDISO
    std::cout << "SOLVING GRID SIZE " << data->hI << " WITH SMOOTHED PARDISO" << std::endl;
    data->SolveWithPreSmoothedPARDISO(1e-10, 0.0, 1000.0, 1);
#else
#ifdef USE_PURE_PARDISO
    std::cout << "SOLVING GRID SIZE " << data->hI << " WITH PARDISO" << std::endl;
    data->SolveWithPARDISO(1e-10, 0.0, 1000.0, 1);
#else
#ifdef USE_PRESMOOTHED_MULTIGRID
    std::cout << "SOLVING GRID SIZE " << data->hI << " WITH PRE-SMOOTHED MULTIGRID" << std::endl;
    data->SolveWithMultigrid(1e-10, 0.0, 1000.0, 1, 8, 1);
#else
#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
    std::cout << "SOLVING GRID SIZE " << data->hI << " WITH PRE-SMOOTHED GAUSS SEIDEL" << std::endl;
    data->SolveWithPseudoMultigrid(1e-10, 0.0, 1000.0, 20000, 1);
#else
#ifdef USE_PURE_GAUSS_SEIDEL
    std::cout << "SOLVING GRID SIZE " << data->hI << " WITH PURE GAUSS SEIDEL" << std::endl;
    data->SolveWithGaussSeidel(1e-10, 0.0, 500.0, 2000, 1);
#endif
#endif
#endif
#endif
#endif
}