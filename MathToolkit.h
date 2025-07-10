#pragma once
#include <cstdint>
namespace MathToolkit
{
	inline double SignedTriangleArea(double x1, double y1, double x2, double y2, double x3, double y3)
	{
		//Area = (1/2) [x1 (y2 – y3) + x2 (y3 – y1) + x3 (y1 – y2)]
		return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
	}
	inline bool ExtractBitFromByte(uint8_t t, int which)
	{
		return ((t >> which) & 1) != 0;
	}
	inline uint32_t CalcRGB(uint8_t r, uint8_t g, uint8_t b)
	{
		uint32_t i = ((0xFF & ((uint32_t) r)) << 16) |
			((0xFF & ((uint32_t) g)) << 8) |
			(0xFF & ((uint32_t) b));
		return i;
	}
	inline uint8_t ExtractRed(uint32_t rgb)
	{
		int red = (rgb >> 16);
		return (uint8_t) red;
	}
	inline uint8_t ExtractGreen(uint32_t rgb)
	{
		int green = (rgb >> 8);
		return (uint8_t) green;
	}
	inline uint8_t ExtractBlue(uint32_t rgb)
	{
		int blue = rgb;
		return (uint8_t) blue;
	}
	inline bool AlmostEquals(int i, double d)
	{
		double diff = d - i;
		return diff <= 0.5 && diff >= -0.5;
	}
	inline int BinaryLog(uint32_t num)
	{
		int log = 0;
		if ((num & 0xffff000) != 0) { num>>=16; log = 16; }
		if (num >= 256) { num >>= 8; log += 8; }
		if (num >= 16) { num >>= 4; log += 4; }
		if (num >= 4) { num >>= 2; log += 2; }
		return log + (num >> 1);
	}
	inline long FastExp(long base, long exponent)
	{
		long answer = 1;
		long powerOfTwo = base;
		for (int i = 0; i < 64; i++)
		{
			if ((exponent & 1) == 1)
				answer *= powerOfTwo;
			powerOfTwo *= powerOfTwo;
			exponent = exponent >> 1;
		}
		return answer;
	}
	inline uint32_t SmoothColorLerp(uint32_t a, uint32_t b, double t)
	{
		if (t <= 0)
			return a;
		if (t >= 1)
			return b;
		uint32_t aR = (a >> 16) & 0xFF, aG = (a >> 8) & 0xFF, aB = a & 0xFF;
		uint32_t bR = (b >> 16) & 0xFF, bG = (b >> 8) & 0xFF, bB = b & 0xFF;

		double red = (1.0 * bR - 1.0 * aR) * (3 - t * 2) * t * t + 1.0 * aR;
		double gre = (1.0 * bG - 1.0 * aG) * (3 - t * 2) * t * t + 1.0 * aG;
		double blu = (1.0 * bB - 1.0 * aB) * (3 - t * 2) * t * t + 1.0 * aB;

		uint32_t R = (uint32_t) red;
		uint32_t G = (uint32_t) gre;
		uint32_t B = (uint32_t) blu;
		return (R << 16) | (G << 8) | B;
	}
	inline double SmoothLerp(double a, double b, double t)
	{
		if (t <= 0)
			return a;
		if (t >= 1)
			return b;
		return (b - a) * (3 - t * 2) * t * t + a;
	}
	inline float SmoothLerp(float a, float b, double t)
	{
		if (t <= 0)
			return a;
		if (t >= 1)
			return b;
		return (float)((b - a) * (3 - t * 2) * t * t + a);
	}
	inline double DerivativeSmoothLerp(double a, double b, double t, double dadt, double dbdt)
	{
		//f(t) = t^2(3-2t)(b(t)-a(t))+a(t)
		//f(t) = (3t^2-2t^3)(b(t)-a(t))+a(t)
		//f'(t) = (3t^2-2t^3)(dbdt-dadt)+(6t-6t^2)(b(t)-a(t))+dadt
		if (t < 0)
			return 0;
		if (t > 1)
			return 0;
		if (t == 0)
			return dadt;
		if (t == 1)
			return dbdt;
		return (3 * t * t - 2 * t * t * t) * (dbdt - dadt) + (6 * t - 6 * t * t) * (b - a) + dadt;
	}
}