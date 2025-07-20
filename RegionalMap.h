#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <glm.hpp>
#include "WatermapValue.h"
#include "LocalMap.h"

class WorldMap;
class SamplePoint;

/*
 * A RegionalMap represents an area 200km x 200km.
 * It is composed of 20x20 squares, each of which is 10km
 * x 10km large.
 */
class RegionalMap
{
public:
	RegionalMap(int x, int y, WorldMap* p);
	~RegionalMap();

	SamplePoint* GetSamplePoint(int index);
	void CalculateAllVoronoiAdjacencies();

	bool IsRidge(double wX, double wY);
	double DrainageBasedRiverPercent(double wX, double wY);
	double RiverPercent(double wX, double wY);

	std::unordered_map<LocalMap*, bool> PrepareForExtensiveEditing();
	void RunFullPhasedErosion();
	void RunFullLaplacianErosion();
	void RunFullRain();

	WatermapValue IsWaterPoint(double wX, double wY);
	double CalculateBaseSedimentDepth(double wX, double wY);
	double CalculateElevation(double wX, double wY);

	double GetElevation(double wX, double wY);
	double GetElevationLaplacian(double wX, double wY);
	glm::dvec2 GetElevationGradient(double wX, double wY);

	LocalMap* GetLocalMapAt(int x, int y);
	LocalMap::Coordinate GetLocalMapAt(double x, double y);
	RegionalMap* GetNorth();
	RegionalMap* GetWest();
	RegionalMap* GetSouth();
	RegionalMap* GetEast();
	int GetOriginX() { return x; }
	int GetOriginY() { return y; }
	int GetWorldX() { return x + ORIGIN_OFFSET; }
	int GetWorldY() { return y + ORIGIN_OFFSET; }
	void EnableRendering() { readyToRender = true; }

	void GatherLocalMapOutlineLocations(
		int tileD,
		int screenWidth, int screenHeight,
		float myX, float myY,
		std::vector<glm::vec2>& outlines);
	void GatherSamplePointCenters(
		float regionalDim,
		int screenWidth, int screenHeight,
		float myX, float myY,
		std::vector<glm::vec2>& pointLocs);

	const std::vector<SamplePoint*>& GetAllPoints() { return voronoiList; }
	SamplePoint* GetAt(int i, int j) { return terrainCells[i * VORONOI_DIM + j]; }
	bool SetPoint(SamplePoint* p);

	std::vector<SamplePoint*> GetNearestN(double wX, double wY, int n);
	SamplePoint* GetNearest(double wX, double wY);
	bool ExistingPointNearby(SamplePoint* cand);

	struct Coordinate
	{
		Coordinate(double x, double y, RegionalMap* owner) : x(x), y(y), src(owner) {}
		double x, y;
		RegionalMap* src;
	};

private:
	double RiverFlowPercent(double flow, double distanceFromFlow);

	double SedimentStairCalculation(
		double terrainHeight, 
		double stepHeight, 
		double flatLength, 
		double flatness, 
		int numSubsteps);

	SamplePoint* GetCell(int i, int j);

public:

	static void ScrollOriginOffsetForOptimalCoastliness();

	inline static constexpr int VORONOI_DIM = 240;
	inline static constexpr double MIN_VORONOI_DIST = 1.41421356237 / VORONOI_DIM;
	inline static constexpr int DIMENSION = 20;
	//if a RegionalMap gets negative coordinates, bad things happen to Perlin RNG
	inline static int ORIGIN_OFFSET = 500;
	inline static constexpr const char* K_HEIGHTMAP_FOLDER_NAME = "Local_Heights_";
	inline static constexpr const char* K_WATERMAP_FOLDER_NAME = "Local_Watermaps_";
	inline static constexpr const char* K_RAINFLOWMAP_FOLDER_NAME = "Local_Rainflowmaps_";
	inline static constexpr const char* K_SEDIMENTMAP_FOLDER_NAME = "Local_Sedimentmaps_";

private:
	int x, y; //while tile am I in the World Map?
	LocalMap* topography[DIMENSION * DIMENSION] = {};
	SamplePoint* terrainCells[VORONOI_DIM * VORONOI_DIM] = {};
	std::vector<SamplePoint*> voronoiList;
	WorldMap* parent;
	bool readyToRender;
};