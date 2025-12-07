#pragma once

#include <glm.hpp>
#include "TileSizes.h"

class RegionalMap;

class LocalMap
{

public:
	LocalMap(int x, int y, RegionalMap* parent);

	bool PrepareForEditing(bool height, bool water, bool rainflow);
	void CompleteEditing(bool height, bool water, bool rainflow, bool killWaterHeightsRec);
	void LaplacianErosionIteration(int num);
	void SendRandomRainErosion(int numDroplets);

	double GetHeight(int px, int py);
	double GetHeight(double lX, double lY);
	glm::dvec2 GetHeightGradient(double lX, double lY);
	double GetHeightLaplacian(double lX, double lY);

	int GetWorldX();
	int GetWorldY();
	LocalMap* GetNorth();
	LocalMap* GetEast();
	LocalMap* GetSouth();
	LocalMap* GetWest();
	LocalMap* GetNorthWest();
	LocalMap* GetNorthEast();
	LocalMap* GetSouthWest();
	LocalMap* GetSouthEast();

public:
	struct Pixel
	{

	};
	struct Coordinate
	{
		Coordinate(double x, double y, LocalMap* owner) : x(x), y(y), src(owner) {}
		double x, y;
		LocalMap* src;
	};

private:
	int x;
	int y;
	RegionalMap* parent;
	bool activityFlag;
public:
	inline static constexpr int METER_DIM = TileSizes::METERS_PER_LOCAL_MAP;
};

