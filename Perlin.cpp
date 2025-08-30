#include "Perlin.h"
#include <cmath>
#include <numbers>
#include "MathToolkit.h"

using namespace std;

//------------------------------------------------------------------------------

double 
Perlin::GetPercentAboveThreshold(double wX, double wY, int octaves) const
{
	double val = Get(wX, wY, octaves);
	val -= threshold;
	val /= (1 - threshold);
	if (val > 1)
		val = 1;
	if (val < 0)
		val = 0;
	return val;
}

//------------------------------------------------------------------------------

double 
Perlin::GetPercentBeneathThreshold(double wX, double wY, int octaves) const
{
	double val = Get(wX, wY, octaves);
	val -= threshold;
	val /= (1 + threshold);
	if (val > 1)
		val = 1;
	if (val < 0)
		val = 0;
	return val;
}

//------------------------------------------------------------------------------

double 
Perlin::GetPercentAboveThreshold(double wX, double wY) const
{
	double val = Get(wX, wY);
	if (absThreshold)
		val = abs(val);
	val -= threshold;
	val /= (1 - threshold);
	if (val > 1)
		val = 1;
	if (val < 0)
		val = 0;
	return val;
}

//------------------------------------------------------------------------------

double
Perlin::GetPercentBeneathThreshold(double wX, double wY) const
{
	double val = Get(wX, wY);
	if (absThreshold)
	{
		val = abs(val);
		val /= (threshold);
		val = 1 - val;
		if (val > 1)
			val = 1;
		if (val < 0)
			val = 0;
		return val;
	}
	else
	{
		val -= threshold;
		val /= (1 + threshold);
		val *= -1;
		if (val > 1)
			val = 1;
		if (val < 0)
			val = 0;
		return val;
	}

}

//------------------------------------------------------------------------------

bool 
Perlin::UnderThreshold(double wX, double wY) const
{
	if (!thresholdSet)
		return true;

	double modThreshold = threshold * GetMaxValue() / NORM_SCALE;
	wX *= scale;
	wY *= scale;
	double result = 0;
	double cM = 1;
	double divider = 1;
	double maxDelta = GetMaxValue();
	for (int i = 0; i < octaves; i++)
	{
		result += PerlinCalc(wX * cM, wY * cM) * divider;
		cM *= 2;
		divider /= octaveScale;
		maxDelta -= divider;

		if (i + 1 < octaves)
		{
			if (absThreshold)
			{
				if (abs(result) + maxDelta < modThreshold)
				{
					return true;
				}
			}
			else
			{
				if (result + maxDelta < modThreshold)
				{
					return true;
				}
			}
		}
	}
	if (absThreshold)
		result = abs(result);

	return result < modThreshold;
}

//------------------------------------------------------------------------------

double 
Perlin::Get(double wX, double wY, int octaves) const
{
	double res = GetNoise(offset + wX * scale, offset + wY * scale, octaves, octaveScale);
	res *= NORM_SCALE;
	if (res >= 1)
	{
		res = 0.99999999;
		//System.out.println("Thought you should know we needed to clamp");
	}
	if (res < -1)
	{
		res = -0.99999999;
		//System.out.println("Thought you should know we needed to clamp");
	}
	return res;
}

//------------------------------------------------------------------------------

double
Perlin::Get(double wX, double wY) const
{
	return Get(wX, wY, octaves);
}

//------------------------------------------------------------------------------

glm::dvec2 
Perlin::GetGradient(double wX, double wY) const
{
	glm::dvec2 res = GetNoiseGradient(offset + wX * scale, offset + wY * scale, octaves, octaveScale);
	res *= scale;
	res *= NORM_SCALE;
	return res;
}

//------------------------------------------------------------------------------

/*private double GetMaxValue(int octave)
{
	double r = 1 / octaveScale;
	double exp = Math.pow(r, octave);
	double n = 1 - exp;
	double d = 1 - r;
	return n / d;
}*/

//------------------------------------------------------------------------------

double
Perlin::GetMaxValue() const
{
	double r = 1 / octaveScale;
	return 1 / (1 - r);
}

//------------------------------------------------------------------------------

double 
Perlin::GetNoise(double x, double y, int octaves, double scaling) const
{
	double result = 0;
	double cM = 1;
	double divider = 1;
	for (int i = 0; i < octaves; i++)
	{
		result += PerlinCalc(x * cM, y * cM) * divider;
		cM *= 2;
		divider /= scaling;
	}
	return result / GetMaxValue();
}

//------------------------------------------------------------------------------

glm::dvec2 
Perlin::GetNoiseGradient(double x, double y, int octaves, double scaling) const
{
	glm::dvec2 result(0, 0);
	double cM = 1;
	double divider = 1;
	for (int i = 0; i < octaves; i++)
	{
		glm::dvec2 octGrad = PerlinGradient(x * cM, y * cM);
		octGrad *= cM;
		octGrad *= divider;
		result += octGrad;
		cM *= 2;
		divider /= scaling;
	}
	result /= GetMaxValue();
	return result;
}

//------------------------------------------------------------------------------

glm::dvec2
Perlin::PerlinGradient(double x, double y) const
{
	int x0 = (int)x;
	int x1 = x0 + 1;
	int y0 = (int)y;
	int y1 = y0 + 1;
	double sX = x - x0;
	double sY = y - y0;

	double n0, n1, n2, n3;
	glm::dvec2 gradN0, gradN1, gradN2, gradN3;

	n0 = DotGridGradientForPerlin(x0, y0, x, y);
	n1 = DotGridGradientForPerlin(x1, y0, x, y);
	n2 = DotGridGradientForPerlin(x0, y1, x, y);
	n3 = DotGridGradientForPerlin(x1, y1, x, y);
	gradN0 = GradDotGridGradientForPerlin(x0, y0, x, y);
	gradN1 = GradDotGridGradientForPerlin(x1, y0, x, y);
	gradN2 = GradDotGridGradientForPerlin(x0, y1, x, y);
	gradN3 = GradDotGridGradientForPerlin(x1, y1, x, y);

	double dx0 = MathToolkit::DerivativeSmoothLerp(n0, n1, sX, gradN0.x, gradN1.x);
	double dx1 = MathToolkit::DerivativeSmoothLerp(n2, n3, sX, gradN2.x, gradN3.x);
	double dx = MathToolkit::SmoothLerp(dx0, dx1, sY);

	double dy0 = MathToolkit::DerivativeSmoothLerp(n0, n2, sY, gradN0.y, gradN2.y);
	double dy1 = MathToolkit::DerivativeSmoothLerp(n1, n3, sY, gradN1.y, gradN3.y);
	double dy = MathToolkit::SmoothLerp(dy0, dy1, sX);

	dx *= 2;
	dx /= 1.41422;
	dy *= 2;
	dy /= 1.41422;

	return glm::dvec2(dx, dy);
}

//------------------------------------------------------------------------------

double
Perlin::PerlinCalc(double x, double y) const
{
	int x0 = (int)x;
	int x1 = x0 + 1;
	int y0 = (int)y;
	int y1 = y0 + 1;
	double sX = x - x0;
	double sY = y - y0;

	double n0, n1, n2, n3, t0, t1;

	n0 = DotGridGradientForPerlin(x0, y0, x, y);
	n1 = DotGridGradientForPerlin(x1, y0, x, y);
	t0 = MathToolkit::SmoothLerp(n0, n1, sX);
	n2 = DotGridGradientForPerlin(x0, y1, x, y);
	n3 = DotGridGradientForPerlin(x1, y1, x, y);
	t1 = MathToolkit::SmoothLerp(n2, n3, sX);

	double result = MathToolkit::SmoothLerp(t0, t1, sY);
	result *= 2;
	result /= 1.41422;
	return result;
}

//------------------------------------------------------------------------------

double 
Perlin::DotGridGradientForPerlin(int ix, int iy, double x, double y) const
{
	double gx = 0, gy = 0;
	uint32_t lx = (uint32_t) ix, ly = (uint32_t) iy;
	//lx += 0x80000000;
	//ly += 0x80000000;
	//lx = ix < 0 ? Integer.MAX_VALUE + ix : ix;
	//ly = iy < 0 ? Integer.MAX_VALUE + iy : iy;
	minstd_rand randOne(lx << 16 | lx >> 16);
	minstd_rand randTwo(ly);
	uint32_t seedX = randOne();
	uint32_t seedY = randTwo();
	uint32_t realSeed = seedX ^ seedY ^ seed;

	minstd_rand newRand(realSeed);
	newRand.discard(5);
	uniform_real_distribution dist(0.0, 1.0);
	double theta = 2 * numbers::pi * dist(newRand);
	gx = (double) sin(theta);
	gy = (double) cos(theta);
	double dx = x - ix;
	double dy = y - iy;
	return gx * dx + gy * dy;
}

//------------------------------------------------------------------------------

glm::dvec2
Perlin::GradDotGridGradientForPerlin(int ix, int iy, double x, double y) const
{
	double gx = 0, gy = 0;
	uint32_t lx = (uint32_t)x, ly = (uint32_t)iy;
	lx += 0x80000000;
	ly += 0x80000000;
	//lx = ix < 0 ? Integer.MAX_VALUE + ix : ix;
	//ly = iy < 0 ? Integer.MAX_VALUE + iy : iy;
	minstd_rand randOne(0xffffffff - lx);
	minstd_rand randTwo(ly);
	uint32_t realSeed = randOne() ^ randTwo() ^ seed;

	minstd_rand newRand(realSeed);
	uniform_real_distribution dist(0.0, 1.0);
	double theta = 2 * numbers::pi * dist(newRand);
	gx = (double) sin(theta);
	gy = (double) cos(theta);
	return glm::dvec2(0, 0);
}

//------------------------------------------------------------------------------

//Determines where land is
const Perlin Perlin::oceans(0.88, 0.0, 10, 1.8, -0.3, false); //0.144 for continents, 0.88 for islands

//Draws strips of mountains
const Perlin Perlin::mountains(0.24, 0, 8, 1.8, 0.036, true); //previously: 0.052

//Marks the center of mountain strips as high peaks
const Perlin Perlin::peaks(Perlin::mountains, 0.007, true); //previously: 0.01

//Creates a ring of hills around mountain strips
const Perlin Perlin::foothills(Perlin::mountains, 0.073, true); //previously: 0.09

//Draws blobs of hills
const Perlin Perlin::randomHills(7.52, 0, 7, 2, 0.18, false);

//Draws random lakes; increasing the threshold decreases the size
const Perlin Perlin::randomLakes(0.88, 0, 7, 2, 0.46, false);

//Cuts through mountain strips, demoting peaks to mountains, mountains to foothills, and foothills to flatland
const Perlin Perlin::randomPasses(1.76, 0, 6, 2, 0.03, true);

//Determines how aggressively we climb
const Perlin Perlin::minMaxSelector(27.04, 0, 4, 1.8, 0, false);

//Adjusts Tectonic Uplift Values
const Perlin Perlin::upliftAdjust(3.6, 0, 10, 2, 0, false);

//Pushes Terrain around a bit to destroy crisp lines
const Perlin Perlin::blurX(80, 0, 4, 2, 0, false);
const Perlin Perlin::blurY(80, 0, 4, 2, 0, false);

const Perlin Perlin::terra_incognita(0.16, 0, 3, 2);

//During erosion simulation, determines how easily the ground erodes
const Perlin Perlin::sedimentStepDelta(34.4, 0, 4, 1.8, 0, false);
const Perlin Perlin::sedimentStepMask(45.6, 0, 6, 1.8, 0.23, true);

//some functions to push the terrain a bit more
const Perlin Perlin::rockyJitters(125, 0, 4, 2, 0, false);

const Perlin Perlin::mountainHeightDelta(172, 0, 6, 2.2, 0, false);
const Perlin Perlin::plainsHeightDelta(88, 0, 6, 2, 0, false);
const Perlin Perlin::elevDeltas[2] = { Perlin::plainsHeightDelta, Perlin::mountainHeightDelta };

//------------------------------------------------------------------------------

bool 
Perlin::SaveSeeds(FILE* wr)
{
	fprintf(wr, "%u\n", oceans.seed);
	fprintf(wr, "%u\n", mountains.seed);
	fprintf(wr, "%u\n", peaks.seed);
	fprintf(wr, "%u\n", foothills.seed);
	fprintf(wr, "%u\n", randomHills.seed);
	fprintf(wr, "%u\n", randomLakes.seed);
	fprintf(wr, "%u\n", randomPasses.seed);
	fprintf(wr, "%u\n", minMaxSelector.seed);
	fprintf(wr, "%u\n", upliftAdjust.seed);
	fprintf(wr, "%u\n", blurX.seed);
	fprintf(wr, "%u\n", blurY.seed);
	fprintf(wr, "%u\n", terra_incognita.seed);
	fprintf(wr, "%u\n", sedimentStepDelta.seed);
	fprintf(wr, "%u\n", sedimentStepMask.seed);
	fprintf(wr, "%u\n", rockyJitters.seed);
	for (Perlin pf : elevDeltas)
		fprintf(wr, "%u\n", pf.seed);
	if (ferror(wr))
		return false;
	return true;
}

//------------------------------------------------------------------------------

bool Perlin::LoadSeeds(std::ifstream& rd)
{
    uint32_t tempSeed = 0;
    if (!(rd >> tempSeed)) return false; oceans.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; mountains.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; peaks.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; foothills.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; randomHills.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; randomLakes.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; randomPasses.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; minMaxSelector.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; upliftAdjust.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; blurX.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; blurY.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; terra_incognita.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; sedimentStepDelta.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; sedimentStepMask.seed = tempSeed;
    if (!(rd >> tempSeed)) return false; rockyJitters.seed = tempSeed;

    for (const Perlin& pf : elevDeltas) {
        if (!(rd >> tempSeed)) return false;
        pf.seed = tempSeed;
    }

    return rd.good();
}
