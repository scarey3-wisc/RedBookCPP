// HydrologyApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <numbers>

#include "..\HydrologySolver\SolverData.h"

template<int e>
void
Visualize(std::filesystem::path outPath, SolverData<e>& data)
{
    std::ofstream out(outPath.c_str());
    out << "P3\n" << data.w << " " << data.h << "\n255\n";
    double minB = std::numeric_limits<double>::max();
    double maxB = std::numeric_limits<double>::min();
    double minW = std::numeric_limits<double>::max();
    double maxW = std::numeric_limits<double>::min();
    for (int y = 0; y < data.h; ++y) {
        for (int x = 0; x < data.w; ++x) {
            if (data.GetW(x, y) < minW)
                minW = data.GetW(x, y);
            if (data.GetW(x, y) > maxW)
                maxW = data.GetW(x, y);
            if (data.GetB(x, y) < minB)
                minB = data.GetB(x, y);
            if (data.GetB(x, y) > maxB)
                maxB = data.GetB(x, y);
        }
    }
    minW = 0;
    std::cout << "B range: " << minB << " - " << maxB << ", W range: " << minW << " - " << maxW << std::endl;
    for (int y = 0; y < data.h; ++y) {
        for (int x = 0; x < data.w; ++x) {
            double w = (data.GetW(x, y) - minW) / (maxW - minW);
            double t = 0;
            if(maxB > minB)
                t = (data.GetB(x, y) - minB) / (maxB - minB);
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

int main(int argc, char** argv)
{

    std::filesystem::path exePath = std::filesystem::absolute(argv[0]);
    std::filesystem::path exeDir = exePath.parent_path();
    std::filesystem::path topDir = exeDir.parent_path().parent_path().parent_path();
    std::filesystem::path outPath = topDir;
    outPath.append("test.ppm");
    outPath.replace_filename("results.ppm");
    SolverData<6> myData(9.8, 0.1, 40);

    for (int j = -1; j <= myData.h; j++)
    {
        for (int i = -1; i <= myData.w; i++)
        {
            double base = 1.0;
            double mainValley = 2.0 * std::abs(myData.w * 0.5 - i) / myData.w;
            double slightSlope = 0.7 * (j + 1.0) / myData.h;
            double sideValley = 0.3 * std::sin(12 * std::numbers::pi * j / myData.h);
            sideValley *= mainValley;
            myData.SetB(i, j, base + mainValley + slightSlope + sideValley);
            myData.SetW(i, j, 0);
        }
    }
    for (int j = 0; j < myData.h; j++)
    {
        for (int i = 0; i < myData.w; i++)
        {
            myData.SetR(i, j, 0.000000049);
            myData.SetW(i, j, 1);
        }
    }

    myData.SolveWithGaussSeidel(1e-10, 0, 10000, 1);
    Visualize<6>(outPath, myData);
}
