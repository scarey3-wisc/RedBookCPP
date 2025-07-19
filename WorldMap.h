#pragma once
#include <string>
#include <vector>
#include "RegionalMap.h"
#include <glad/glad.h>
#include "OutlineRendering.h"


class SamplePoint;

class WorldMap
{

public:

	WorldMap(const char* name);
	~WorldMap();

	void EnableAllRegionalRenderings();
	RegionalMap* GetRegion(int x, int y);
	RegionalMap* RequestRegion(int x, int y);
	bool CreateRegion(int x, int y);
	void ExpandEmptySpace(int posXInc, int negXInc, int posYInc, int negYInc);

	RegionalMap::Coordinate GetRegionalMapAt(double x, double y);
	RegionalMap::Coordinate GetRegionalMapAt(double x, double y, bool createIfNot);

	void InitPoissonDiscSample(RegionalMap* target);
	std::vector<SamplePoint*> PoissonDiscSample(std::vector<SamplePoint*>& active);
	void FillAllContinents(int regionX, int regionY, std::vector<SamplePoint*>& allFound);

	double GetTileSize() { return tileSize; }

	bool TextureInit() const { return textureInit; }
	GLuint Render(int width, int height);

private:
	void InitNewWorld();

	void InitializeRenderTarget(int width, int height);

public:
	bool FileAvaiable;

	double dX, dY;
	int tileSize;

private:
	int x0, y0; //which cell in regions counts as "0,0"? This may change as we expand!
	int w, h; //how many cells do we have?
	RegionalMap** regions;
	//WorldMapTool activeTool;
	//ZoomUpdater myZoom;
	std::string worldName;


private:
	GLuint frameBuffer = 0;
	GLuint renderTexture = 0;
	OutlineRendering outlineRenderer;

	bool textureInit = false;
	int texWidth;
	int texHeight;

public:
	inline static constexpr int DEFAULT_TILE_SIZE = 10; //how many pixels does a local map get?
	inline static constexpr int MINIMUM_TILE_SIZE = 9;
	inline static constexpr int MAXIMUM_TILE_SIZE = 2047; //1023;
	inline static constexpr const char* K_EMPTY_REGION_NAME = "null_region";
	inline static constexpr const char* K_SAVE_FOLDER_NAME = "Worlds";
	inline static constexpr const char* K_REGIONS_FOLDER_NAME = "Regions";
	inline static constexpr int SCALE_LABEL_TARGET_WIDTH = 200;
	inline static constexpr const char* TARGET_SCALE_LABLES[] =
	{
		"200 Miles",
		"100 Miles",
		"50 Miles",
		"25 Miles",
		"10 Miles",
		"5 Miles",
		"2 Miles",
		"1 Mile",
		"1/2 Mile",
		"1/3 Mile",
		"1/4 Mile",
		"1000 Feet",
		"500 Feet",
		"200 Feet",
		"100 Feet",
		"50 Feet",
		"20 Feet",
		"10 Feet"
	};
	inline static constexpr double LABEL_DISTANCES[] =
	{
		321869,
		160934,
		80467.2,
		40233.6,
		16093.4,
		8046.72,
		3218.69,
		1609.34,
		804.672,
		536.448,
		406.336,
		304.8,
		152.4,
		60.96,
		30.48,
		15.24,
		6.096,
		3.048
	};
};

