// HydrologyApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "..\HydrologySolver\Globals.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <numbers>

#include "..\HydrologySolver\SolverData.h"

template<int e>
void
Visualize(std::filesystem::path outPath, SolverData<e>& data, int mode)
{
    std::ofstream out(outPath.c_str());
    out << "P3\n" << data.wO << " " << data.hO << "\n255\n";
    double minB = std::numeric_limits<double>::max();
    double maxB = std::numeric_limits<double>::min();
    double minW = std::numeric_limits<double>::max();
    double maxW = std::numeric_limits<double>::min();
    for (int y = 0; y < data.hO; ++y) {
        for (int x = 0; x < data.wO; ++x) {
            if (data.GetW(x, y) < minW)
                minW = data.GetW(x, y);
            if (data.GetW(x, y) > maxW)
                maxW = data.GetW(x, y);
            if (mode == 2)
            {
                if (data.GetZ(x, y) < minB)
                    minB = data.GetZ(x, y);
                if (data.GetZ(x, y) > maxB)
                    maxB = data.GetZ(x, y);
            }
            else
            {
                if (data.GetB(x, y) < minB)
                    minB = data.GetB(x, y);
                if (data.GetB(x, y) > maxB)
                    maxB = data.GetB(x, y);
            }

        }
    }
    minW = 0;
    std::cout << "B range: " << minB << " - " << maxB << ", W range: " << minW << " - " << maxW << std::endl;
    for (int y = 0; y < data.hO; ++y) {
        for (int x = 0; x < data.wO; ++x) {
            double w = (data.GetW(x, y) - minW) / (maxW - minW);
            double t = 0;
            if(maxB > minB)
                t = (data.GetB(x, y) - minB) / (maxB - minB);
            if(mode == 2 && maxB > minB)
				t = (data.GetZ(x, y) - minB) / (maxB - minB);
            if (mode == 0)
                w = std::pow(w, 0.5);
            else
                w = 0;
            /*
            int wR = 19;
            int wG = 113;
            int wB = 176;
            int lR = 105;
            int lG = 73;
            int lB = 35;
            int hR = 209;
            int hG = 199;
            int hB = 174;
            */
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

static constexpr int dim = 10;

int main(int argc, char** argv)
{
    std::filesystem::path exePath = std::filesystem::absolute(argv[0]);
    std::filesystem::path exeDir = exePath.parent_path();
    std::filesystem::path topDir = exeDir.parent_path().parent_path().parent_path();
    std::filesystem::path outPath = topDir;
    MY_PATH = outPath;
	MY_PATH.append("output");

    SolverData<dim> myData(9.8, 0.1, 40 * 1023);

    for (int j = 0; j < myData.hO; j++)
    {
        for (int i = 0; i < myData.wO; i++)
        {
            double base = 0.01;
            double mainValley = 40.0 * std::abs(myData.wO * 0.5 - i) / myData.wO;
            double slightSlope = 40.0 * (j + 1.0) / myData.hO;
            double sideValley = 0.5 + 0.5 * std::sin(8 * std::numbers::pi * j / myData.hO);
            sideValley *= 0.07 * mainValley * mainValley;

            double lake = 0;
            if((i > .25 * myData.wO && i < .32 * myData.wO) &&
               (j > .34 * myData.hO && j < .55 * myData.hO))
            {
                lake = -4.0;
			}
            myData.SetB(i, j, base + mainValley + slightSlope + sideValley + lake);
            //myData.SetB(i, j, 1);
            myData.SetW(i, j, 0);
        }
    }
    for (int j = 1; j < myData.hO - 1; j++)
    {
        for (int i = 1; i < myData.wO - 1; i++)
        {
            myData.SetR(i, j, 0.000000049);
            myData.SetW(i, j, 1);
        }
    }
    outPath.append("test.ppm");

#ifdef USE_PARDISO
	std::cout << "SOLVING GRID SIZE " << myData.hI << " WITH PARDISO" << std::endl;
    myData.SolveWithPARDISO(1e-10, 0, 10000.0, 1);
#else
#ifdef USE_PRESMOOTHED_MULTIGRID
    std::cout << "SOLVING GRID SIZE " << myData.hI << " WITH PRE-SMOOTHED MULTIGRID" << std::endl;
    myData.SolveWithMultigrid(1e-10, 0, 10000.0, 2, 2000, 1);
#else
#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
    std::cout << "SOLVING GRID SIZE " << myData.hI << " WITH PRE-SMOOTHED GAUSS SEIDEL" << std::endl;
    myData.SolveWithPseudoMultigrid(1e-10, 0, 10000, 2000, 1);
#else
#ifdef USE_PURE_GAUSS_SEIDEL
    std::cout << "SOLVING GRID SIZE " << myData.hI << " WITH PURE GAUSS SEIDEL" << std::endl;
    myData.SolveWithGaussSeidel(1e-10, 0, 10000.0, 2000, 1);
#endif
#endif
#endif
#endif
	
	
    outPath.replace_filename("waterways.ppm");
    Visualize<dim>(outPath, myData, 0);
    outPath.replace_filename("terrain.ppm");
    Visualize<dim>(outPath, myData, 1);
    outPath.replace_filename("surface.ppm");
    Visualize<dim>(outPath, myData, 2);
}
