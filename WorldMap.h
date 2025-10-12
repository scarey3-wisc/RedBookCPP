#pragma once
#include <string>
#include <vector>
#include "RegionalMap.h"
#include <glad/glad.h>
#include "OutlineRendering.h"
#include "MeshPointRendering.h"
#include "VoronoiRendering.h"
#include "RegionalDataRendering.h"
#include "Heightmap.h"


class SamplePoint;

class WorldMap
{

public:

	WorldMap(const char* name, const std::array<int, REGIONAL_MAP_NUM_LODS>& cacheSizes);
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

	template<typename Func>
	void ForEachOnSceenRegionalMap(Func f, int width, int height)
	{
		float regionDim = (float) (tileSize * RegionalMap::DIMENSION);
		float xOrigin = 0.0f, yOrigin = 0.0f;
		xOrigin += width / 2, yOrigin += height / 2;
		xOrigin -= 1 * x0 * regionDim; yOrigin -= 1 * y0 * regionDim;
		xOrigin -= 0.5f * regionDim, yOrigin -= 0.5f * regionDim;
		xOrigin += (float)dX, yOrigin += (float)dY;
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				bool onscreen = true;
				if (xOrigin > width)
					onscreen = false;
				if (yOrigin > height)
					onscreen = false;
				if (xOrigin + regionDim < 0)
					onscreen = false;
				if (yOrigin + regionDim < 0)
					onscreen = false;
				if (onscreen && regions[i * h + j] != nullptr)
					f(regions[i * h + j], xOrigin, yOrigin);
				yOrigin += regionDim;
			}
			xOrigin += regionDim;
			yOrigin -= h * regionDim;
		}
	}

	void RenderRegionalData(int width, int height, float regionDim);
	void RenderVoronoiTiles(int width, int height, float regionDim, bool terrainBased);
	void RenderOutlines(int width, int height, float regionDim);
	void RenderMeshPoints(int width, int height, float regionDim);

	void InitNewWorld();

	void InitializeRenderTarget(int width, int height);

public:
	bool FileAvaiable;

	double dX, dY;
	int tileSize;
	HeightmapManager heightmaps;

	RegionalDataRequestStack myWorkerThread = RegionalDataRequestStack(1);

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
	GLuint depthBuffer = 0;
	OutlineRendering outlineRenderer;
	MeshPointRendering meshPointRenderer;
	VoronoiRendering voronoiRenderer;
	RegionalDataRendering regionalDataRenderer;

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
};

