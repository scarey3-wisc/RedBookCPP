#pragma once

#include <random>
#include <cstdint>

class Perlin
{
public:
	inline Perlin(Perlin pf, double t, bool abs)
		: scale(pf.scale), offset(pf.offset), octaves(pf.octaves), octaveScale(pf.octaveScale), threshold(t),
		thresholdSet(true), absThreshold(abs), seed(pf.seed)
	{

	}
	inline Perlin(double s, double o, int oct, double fs)
		: scale(s), offset(o), octaves(oct), octaveScale(fs), threshold(0),
		thresholdSet(false), absThreshold(false), seed(0)
	{
		GenerateSeed();
	}
	inline Perlin(double s, double o, int oct, double fs, double t, bool abs)
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

	inline static constexpr double BLUR_DISTANCE = 450;
	inline static constexpr double rockJitterScale = 450;
	inline static constexpr double elevDeltaScales[2] = { 40, 90 };



private:

	double Get(double wX, double wY, int octaves);

	inline void GenerateSeed()
	{
		std::random_device rd;
		std::mt19937_64 eng(rd());
		seed = eng();
	}
	double scale;
	double offset;
	double octaveScale;
	int octaves;
	double threshold;
	bool thresholdSet;
	bool absThreshold;
	uint64_t seed;
	inline static constexpr double NORM_SCALE = 1.5;

};


//Determines where land is
inline static const Perlin perlin_oceans( 1.1, 0.0, 10, 1.8, -0.3, false ); //0.18 for continents, 1.1 for islands
//Draws strips of mountains
inline static const Perlin perlin_mountains(0.3, 0, 8, 1.8, 0.036, true); //previously: 0.052

//Marks the center of mountain strips as high peaks
inline static const Perlin perlin_peaks(perlin_mountains, 0.007, true); //previously: 0.01

//Creates a ring of hills around mountain strips
inline static const Perlin perlin_foothills(perlin_mountains, 0.073, true); //previously: 0.09

//Draws blobs of hills
inline static const Perlin perlin_randomHills(9.4, 0, 7, 2, 0.18, false);

//Draws random lakes; increasing the threshold decreases the size
inline static const Perlin perlin_randomLakes(1.1, 0, 7, 2, 0.46, false);

//Cuts through mountain strips, demoting peaks to mountains, mountains to foothills, and foothills to flatland
inline static const Perlin perlin_randomPasses(2.2, 0, 6, 2, 0.03, true);

//Determines how aggressively we climb
inline static const Perlin perlin_minMaxSelector(33.8, 0, 4, 1.8, 0, false);

//Adjusts Tectonic Uplift Values
inline static const Perlin perlin_upliftAdjust(4.5, 0, 10, 2, 0, false);

//Pushes Terrain around a bit to destroy crisp lines
inline static const Perlin perlin_blurX(100, 0, 4, 2, 0, false);
inline static const Perlin perlin_blurY(100, 0, 4, 2, 0, false);

inline static const Perlin perlin_terra_incognita(0.1, 0, 3, 2);

//During erosion simulation, determines how easily the ground erodes
inline static const Perlin perlin_sedimentStepDelta(43, 0, 4, 1.8, 0, false);
inline static const Perlin perlin_sedimentStepMask(57, 0, 6, 1.8, 0.23, true);

//some functions to push the terrain a bit more
inline static const Perlin perlin_rockyJitters(250, 0, 4, 2, 0, false);

inline static const Perlin perlin_mountainHeightDelta(215, 0, 6, 2.2, 0, false);
inline static const Perlin perlin_plainsHeightDelta(110, 0, 6, 2, 0, false);
inline static const Perlin perlin_elevDeltas[2] = { perlin_plainsHeightDelta, perlin_mountainHeightDelta };