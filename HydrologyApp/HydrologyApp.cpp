// HydrologyApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "..\HydrologySolver\FileGlobals.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <numbers>

#include "..\HydrologySolver\SolverData.h"
#include "Timer.h"

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

template<int a, int b>
double 
ErrorNorm(SolverData<a>& coarseSolve, SolverData<b>& fineSolve)
{
    int powDif = b - a;
    int mult = (1 << powDif);
    double norm = 0;
    for (int j = 0; j < coarseSolve.hO; j++)
    {
        for (int i = 0; i < coarseSolve.wO; i++)
        {
            double coarseW = coarseSolve.GetW(i, j);
            double fineW = fineSolve.GetW(i * mult, j * mult);
            double del = fineW - coarseW;
            norm += del * del;
        }
    }
    return sqrt(norm / (coarseSolve.wI * coarseSolve.hI));
}

template<int e>
void SetRain(SolverData<e>& data, double rainScale)
{
    for (int j = 0; j < data.hO; j++)
    {
        for (int i = 0; i < data.wO; i++)
        {
            //this 4.9 * 10^-8 is the annual rainfall New Orleans
            //gets in m / s (because it's something like 64 inches / year)
            data.SetR(i, j, 0.000000049 * rainScale);
        }
    }
}

template<int e>
void SetInitialWaterGuess(SolverData<e>& data)
{
    for (int j = 0; j < data.hO; j++)
    {
        for (int i = 0; i < data.wO; i++)
        {
            data.SetW(i, j, 0);
        }
    }
    for (int j = 1; j < data.hO - 1; j++)
    {
        for (int i = 1; i < data.wO - 1; i++)
        {
            data.SetW(i, j, 1);
        }
    }
}

template<int e>
int LaunchSolve(SolverData<e>& data, bool verbose)
{
    int verbosity = 0;
    if (verbose)
        verbosity = 1;
#ifdef USE_PRESMOOTHED_PARDISO
    if(verbose)
        std::cout << "SOLVING GRID SIZE " << data.hI << " WITH SMOOTHED PARDISO" << std::endl;
    return data.SolveWithPreSmoothedPARDISO(1e-10, 0.0, 10000.0, verbosity);
#else
#ifdef USE_PURE_PARDISO
    if(verbose)
        std::cout << "SOLVING GRID SIZE " << data.hI << " WITH PARDISO" << std::endl;
    return data.SolveWithPARDISO(1e-10, 0.0, 10000.0, verbosity);
#else
#ifdef USE_PRESMOOTHED_MULTIGRID
    if(verbose)
        std::cout << "SOLVING GRID SIZE " << data.hI << " WITH PRE-SMOOTHED MULTIGRID" << std::endl;
    return data.SolveWithMultigrid(1e-10, 0.0, 10000.0, 1, 8, verbosity);
#else
#ifdef USE_PRESMOOTHED_GAUSS_SEIDEL
    if(verbose)
        std::cout << "SOLVING GRID SIZE " << data.hI << " WITH PRE-SMOOTHED GAUSS SEIDEL" << std::endl;
    return data.SolveWithPseudoMultigrid(1e-10, 0.0, 10000.0, 2000, verbosity);
#else
#ifdef USE_PURE_GAUSS_SEIDEL
    if(verbose)
        std::cout << "SOLVING GRID SIZE " << data.hI << " WITH PURE GAUSS SEIDEL" << std::endl;
    return data.SolveWithGaussSeidel(1e-10, 0.0, 10000.0, 2000, verbosity);
#endif
#endif
#endif
#endif
#endif
    return 0;
}

template <int e>
double TimeSolve(SolverData<e>& data, int numTests, bool verbose)
{
    Timer watch;
    double totalTime = 0.0;
    for (int i = 0; i < numTests; i++)
    {
        SetInitialWaterGuess(data);
        watch.Start();
        LaunchSolve(data, false);
        if (verbose)
            totalTime += watch.Stop("Solve " + std::to_string(i) + ": ");
        else
        {
            totalTime += watch.Stop();
            std::cout << 'X';
        }
    }
    if (!verbose)
        std::cout << ' ';
    return totalTime / numTests;
}


template<int e>
void SetSinusoidalValley(SolverData<e>& data)
{
    double base = 0.01;
    for (int j = 0; j < data.hO; j++)
    {
        for (int i = 0; i < data.wO; i++)
        {
            double xNorm = 2.0 * (i - (data.wO - 1) / 2) / (data.wO - 1);
            double yNorm = 2.0 * (j - (data.hO - 1) / 2) / (data.hO - 1);
            //just getting all our coordinates on a [-1 - 1] x [-1 - 1] domain
            double mainValley = 20.0 * std::abs(xNorm);
            double slightSlope = 20.0 * (yNorm + 1);
            double sideValley = 0.5 + 0.5 * std::sin(4 * std::numbers::pi * (yNorm + 1));
            sideValley *= 0.07 * mainValley * mainValley;

            double lake = 0;
            if ((xNorm > -0.5 && xNorm < -0.36) &&
                (yNorm > -0.32 && yNorm < 0.1))
            {
                lake = -4.0;
            }
            double val = base + mainValley + slightSlope + sideValley + lake;
            data.SetB(i, j, val);
        }
    }
    SetRain(data, 1);
}

void Scale_Test_SinusoidalValley()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "SINUSOIDAL VALLEY TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetSinusoidalValley(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << TimeSolve(s5, 5, false) << "ms" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetSinusoidalValley(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << TimeSolve(s6, 5, false) << "ms" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetSinusoidalValley(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << TimeSolve(s7, 5, false) << "ms" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetSinusoidalValley(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << TimeSolve(s8, 5, false) << "ms" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetSinusoidalValley(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << TimeSolve(s9, 5, false) << "ms" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetSinusoidalValley(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << TimeSolve(s10, 5, false) << "ms" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetSinusoidalValley(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << TimeSolve(s11, 5, false) << "ms" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetSinusoidalValley(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << TimeSolve(s12, 5, false) << "ms" << std::endl;
}

void Itr_Test_SinusoidalValley()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "SINUSOIDAL VALLEY TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetSinusoidalValley(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << LaunchSolve(s5, false) << " itr" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetSinusoidalValley(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << LaunchSolve(s6, false) << " itr" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetSinusoidalValley(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << LaunchSolve(s7, false) << " itr" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetSinusoidalValley(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << LaunchSolve(s8, false) << " itr" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetSinusoidalValley(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << LaunchSolve(s9, false) << " itr" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetSinusoidalValley(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << LaunchSolve(s10, false) << " itr" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetSinusoidalValley(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << LaunchSolve(s11, false) << " itr" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetSinusoidalValley(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << LaunchSolve(s12, false) << " itr" << std::endl;
}

void Norm_Test_SinusoidalValley()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "SINUSOIDAL VALLEY TEST" << std::endl << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetSinusoidalValley(s12);
    SetInitialWaterGuess(s12);
    LaunchSolve(s12, false);
    std::cout << "Baseline Calculated" << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetSinusoidalValley(s5);
    SetInitialWaterGuess(s5);
    LaunchSolve(s5, false);
    std::cout << "Mesh Level 5: " << ErrorNorm(s5, s12) << " err" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetSinusoidalValley(s6);
    SetInitialWaterGuess(s6);
    LaunchSolve(s6, false);
    std::cout << "Mesh Level 6: " << ErrorNorm(s6, s12) << " err" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetSinusoidalValley(s7);
    SetInitialWaterGuess(s7);
    LaunchSolve(s7, false);
    std::cout << "Mesh Level 7: " << ErrorNorm(s7, s12) << " err" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetSinusoidalValley(s8);
    SetInitialWaterGuess(s8);
    LaunchSolve(s8, false);
    std::cout << "Mesh Level 8: " << ErrorNorm(s8, s12) << " err" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetSinusoidalValley(s9);
    SetInitialWaterGuess(s9);
    LaunchSolve(s9, false);
    std::cout << "Mesh Level 9: " << ErrorNorm(s9, s12) << " err" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetSinusoidalValley(s10);
    SetInitialWaterGuess(s10);
    LaunchSolve(s10, false);
    std::cout << "Mesh Level 10: " << ErrorNorm(s10, s12) << " err" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetSinusoidalValley(s11);
    SetInitialWaterGuess(s11);
    LaunchSolve(s11, false);
    std::cout << "Mesh Level 11: " << ErrorNorm(s11, s12) << " err" << std::endl;
}

template<int e>
void SetCaldera(SolverData<e>& data)
{
    double base = 0.01;
    for (int j = 0; j < data.hO; j++)
    {
        for (int i = 0; i < data.wO; i++)
        {
            double xNorm = 2.0 * (i - (data.wO - 1) / 2) / (data.wO - 1);
            double yNorm = 2.0 * (j - (data.hO - 1) / 2) / (data.hO - 1);
            //just getting all our coordinates on a [-1 - 1] x [-1 - 1] domain
            double radSqr = xNorm * xNorm + yNorm * yNorm;
            radSqr *= 20 * 20;

            double offsetRadSqr = (xNorm - 0.05) * (xNorm - 0.05) + yNorm * yNorm;
            offsetRadSqr *= 20 * 20;

            double basin = 0.5 * radSqr * std::exp(-0.04 * radSqr);
            double offsetBasin = 0.5 * offsetRadSqr * std::exp(-0.01 * radSqr);
            double val = 4 - basin + offsetBasin;

            data.SetB(i, j, val);
        }
    }
    SetRain(data, 1);
}

void Scale_Test_Caldera()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "CALDERA TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetCaldera(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << TimeSolve(s5, 5, false) << "ms" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetCaldera(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << TimeSolve(s6, 5, false) << "ms" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetCaldera(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << TimeSolve(s7, 5, false) << "ms" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetCaldera(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << TimeSolve(s8, 5, false) << "ms" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetCaldera(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << TimeSolve(s9, 5, false) << "ms" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetCaldera(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << TimeSolve(s10, 5, false) << "ms" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetCaldera(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << TimeSolve(s11, 5, false) << "ms" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetCaldera(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << TimeSolve(s12, 5, false) << "ms" << std::endl;
}

void Itr_Test_Caldera()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "CALDERA TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetCaldera(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << LaunchSolve(s5, false) << " itr" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetCaldera(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << LaunchSolve(s6, false) << " itr" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetCaldera(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << LaunchSolve(s7, false) << " itr" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetCaldera(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << LaunchSolve(s8, false) << " itr" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetCaldera(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << LaunchSolve(s9, false) << " itr" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetCaldera(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << LaunchSolve(s10, false) << " itr" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetCaldera(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << LaunchSolve(s11, false) << " itr" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetCaldera(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << LaunchSolve(s12, false) << " itr" << std::endl;
}

template<int e>
void SetFurrowedSlope(SolverData<e>& data)
{
    double base = 0.01;
    for (int j = 0; j < data.hO; j++)
    {
        for (int i = 0; i < data.wO; i++)
        {
            double xNorm = 2.0 * (i - (data.wO - 1) / 2) / (data.wO - 1);
            double yNorm = 2.0 * (j - (data.hO - 1) / 2) / (data.hO - 1);
            //just getting all our coordinates on a [-1 - 1] x [-1 - 1] domain

            double slopeX = 0.6 * 20 * xNorm;
            double slopeY = 0.3 * 20 * yNorm;
            double waveX = cos(20 * xNorm) + cos(20 * xNorm) * cos(20 * xNorm);
            double pitRad = (xNorm * 20 - 15) * (xNorm * 20 - 15) + yNorm * 20 * yNorm * 20;
            double pit = -11 * exp(-0.1 * pitRad);
            double pit2Rad = (xNorm * 20 + 15) * (xNorm * 20 - 15) + yNorm * 20 * yNorm * 20;
            double pit2 = -18 * exp(-0.01 * pit2Rad);

            double val = 20 + slopeX + slopeY + waveX + pit + pit2;

            data.SetB(i, j, val);
        }
    }
    SetRain(data, 1);
}

void Scale_Test_FurrowedSlope()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "FURROWED SLOPE TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetFurrowedSlope(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << TimeSolve(s5, 5, false) << "ms" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetFurrowedSlope(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << TimeSolve(s6, 5, false) << "ms" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetFurrowedSlope(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << TimeSolve(s7, 5, false) << "ms" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetFurrowedSlope(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << TimeSolve(s8, 5, false) << "ms" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetFurrowedSlope(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << TimeSolve(s9, 5, false) << "ms" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetFurrowedSlope(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << TimeSolve(s10, 5, false) << "ms" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetFurrowedSlope(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << TimeSolve(s11, 5, false) << "ms" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetFurrowedSlope(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << TimeSolve(s12, 5, false) << "ms" << std::endl;
}

void Itr_Test_FurrowedSlope()
{
    double domainSize = 40 * 1024;
    double gravity = 9.8;
    double drag = 0.7;
    std::cout << "FURROWED SLOPE TEST" << std::endl << std::endl;

    SolverData<5> s5(gravity, drag, domainSize);
    SetFurrowedSlope(s5);
    SetInitialWaterGuess(s5);
    std::cout << "Mesh Level 5: " << LaunchSolve(s5, false) << " itr" << std::endl;

    SolverData<6> s6(gravity, drag, domainSize);
    SetFurrowedSlope(s6);
    SetInitialWaterGuess(s6);
    std::cout << "Mesh Level 6: " << LaunchSolve(s6, false) << " itr" << std::endl;

    SolverData<7> s7(gravity, drag, domainSize);
    SetFurrowedSlope(s7);
    SetInitialWaterGuess(s7);
    std::cout << "Mesh Level 7: " << LaunchSolve(s7, false) << " itr" << std::endl;

    SolverData<8> s8(gravity, drag, domainSize);
    SetFurrowedSlope(s8);
    SetInitialWaterGuess(s8);
    std::cout << "Mesh Level 8: " << LaunchSolve(s8, false) << " itr" << std::endl;

    SolverData<9> s9(gravity, drag, domainSize);
    SetFurrowedSlope(s9);
    SetInitialWaterGuess(s9);
    std::cout << "Mesh Level 9: " << LaunchSolve(s9, false) << " itr" << std::endl;

    SolverData<10> s10(gravity, drag, domainSize);
    SetFurrowedSlope(s10);
    SetInitialWaterGuess(s10);
    std::cout << "Mesh Level 10: " << LaunchSolve(s10, false) << " itr" << std::endl;

    SolverData<11> s11(gravity, drag, domainSize);
    SetFurrowedSlope(s11);
    SetInitialWaterGuess(s11);
    std::cout << "Mesh Level 11: " << LaunchSolve(s11, false) << " itr" << std::endl;

    SolverData<12> s12(gravity, drag, domainSize);
    SetFurrowedSlope(s12);
    SetInitialWaterGuess(s12);
    std::cout << "Mesh Level 12: " << LaunchSolve(s12, false) << " itr" << std::endl;
}

template<int dim>
void
main_test_appearance()
{
    SolverData<dim> myData(9.8, 0.7, 40 * 1024);
    SetRain(myData, 1);
    SetInitialWaterGuess(myData);

    //SetSinusoidalValley(myData);
    //SetCaldera(myData);
    SetFurrowedSlope(myData);

    LaunchSolve(myData, true);
    std::filesystem::path outPath = MY_PATH;
    outPath.replace_filename("waterways.ppm");
    Visualize<dim>(outPath, myData, 0);
    outPath.replace_filename("terrain.ppm");
    Visualize<dim>(outPath, myData, 1);
    outPath.replace_filename("surface.ppm");
    Visualize<dim>(outPath, myData, 2);
}

int main(int argc, char** argv)
{
    std::filesystem::path exePath = std::filesystem::absolute(argv[0]);
    std::filesystem::path exeDir = exePath.parent_path();
    std::filesystem::path topDir = exeDir.parent_path().parent_path().parent_path();
    std::filesystem::path outPath = topDir;
    MY_PATH = outPath;
	MY_PATH.append("output");

    main_test_appearance<8>();

    //Scale_Test_SinusoidalValley();
    //Scale_Test_Caldera();
    //Scale_Test_FurrowedSlope();

    //Itr_Test_SinusoidalValley();
    //Itr_Test_Caldera();
    //Itr_Test_FurrowedSlope();


    //Norm_Test_SinusoidalValley();
}
