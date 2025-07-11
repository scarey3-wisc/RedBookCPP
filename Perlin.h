#pragma once

#include <random>
#include <cstdint>
#include <glm.hpp>

class Perlin
{
public:
	Perlin(Perlin pf, double t, bool abs)
		: scale(pf.scale), offset(pf.offset), octaves(pf.octaves), octaveScale(pf.octaveScale), threshold(t),
		thresholdSet(true), absThreshold(abs), seed(pf.seed)
	{

	}
	Perlin(double s, double o, int oct, double fs)
		: scale(s), offset(o), octaves(oct), octaveScale(fs), threshold(0),
		thresholdSet(false), absThreshold(false), seed(0)
	{
		GenerateSeed();
	}
	Perlin(double s, double o, int oct, double fs, double t, bool abs)
		: scale(s), offset(o), octaves(oct), octaveScale(fs), threshold(t),
		thresholdSet(true), absThreshold(abs), seed(0)
	{
		GenerateSeed();
	}

	double GetPercentAboveThreshold(double wX, double wY, int octaves);
	double GetPercentBeneathThreshold(double wX, double wY, int octaves);
	double GetPercentAboveThreshold(double wX, double wY);
	double GetPercentBeneathThreshold(double wX, double wY);
	bool UnderThreshold(double wX, double wY);
	double Get(double wX, double wY);
	glm::dvec2 GetGradient(double wX, double wY);

private:

	double Get(double wX, double wY, int octaves);
	double GetMaxValue();
	double GetNoise(double x, double y, int octaves, double scaling);
	glm::dvec2 GetNoiseGradient(double x, double y, int octaves, double scaling);
	glm::dvec2 PerlinGradient(double x, double y);
	double PerlinCalc(double x, double y);
	double DotGridGradientForPerlin(int ix, int iy, double x, double y);
	glm::dvec2 GradDotGridGradientForPerlin(int ix, int iy, double x, double y);


	void GenerateSeed()
	{
		std::random_device rd;
		std::mt19937 eng(rd());
		seed = eng();
	}
	double scale;
	double offset;
	double octaveScale;
	int octaves;
	double threshold;
	bool thresholdSet;
	bool absThreshold;
	uint32_t seed;
	inline static constexpr double NORM_SCALE = 1.5;

public: 

	inline static constexpr double BLUR_DISTANCE = 450;
	inline static constexpr double rockJitterScale = 450;
	inline static constexpr double elevDeltaScales[2] = { 40, 90 };

	static const Perlin oceans;
	static const Perlin mountains;
	static const Perlin peaks;
	static const Perlin foothills;
	static const Perlin randomHills;
	static const Perlin randomLakes;
	static const Perlin randomPasses;
	static const Perlin minMaxSelector;
	static const Perlin upliftAdjust;
	static const Perlin blurX;
	static const Perlin blurY;
	static const Perlin terra_incognita;
	static const Perlin sedimentStepDelta;
	static const Perlin sedimentStepMask;
	static const Perlin rockyJitters;
	static const Perlin mountainHeightDelta;
	static const Perlin plainsHeightDelta;
	static const Perlin elevDeltas[2];

	static bool SaveSeeds(FILE* wr);
	static bool LoadSeeds(FILE* rd);
};