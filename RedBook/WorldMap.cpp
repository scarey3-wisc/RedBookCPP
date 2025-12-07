#include "WorldMap.h"
#include "RegionalMap.h"
#include "SamplePoint.h"
#include "GlobalRand.h"
#include "Switches.h"
#include <cmath>
#include <numbers>
#include <deque>
#include <unordered_set>
#include <iostream>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>

using namespace std;
WorldMap::WorldMap(const char* name, const std::array<int, REGIONAL_MAP_NUM_LODS>& cacheSizes) :
	worldName(name), FileAvaiable(true), regions(nullptr),
	x0(0), y0(0), w(1), h(1), dX(0), dY(0), texWidth(0), texHeight(0),
	heightmaps(cacheSizes, this, &myWorkerThread), depthmaps(cacheSizes, this, &myWorkerThread)
{
	tileSize = DEFAULT_TILE_SIZE;
	InitNewWorld();
	/*File topDir = new File(GetDirectory());
	if (topDir.exists() && topDir.isDirectory())
	{
		FileAvailable = true;
	}
	else
	{
		FileAvailable = false;
		topDir.mkdir();
		InitNewWorld();
	}*/
}

WorldMap::~WorldMap()
{
	if (regions != nullptr)
	{
		for (int i = 0; i < w * h; i++)
		{
			if (regions[i] != nullptr)
				delete regions[i];
		}
		delete[] regions;
		regions = nullptr;
	}
	if (frameBuffer != 0)
		glDeleteFramebuffers(1, &frameBuffer);
	if (renderTexture != 0)
		glDeleteTextures(1, &renderTexture);

	outlineRenderer.Cleanup();
	meshPointRenderer.Cleanup();
	voronoiRenderer.Cleanup();
	regionalDataRenderer.Cleanup();
}
/*public String GetDirectory()
{
	return K_SAVE_FOLDER_NAME + File.separator + worldName;
}
public void DeleteWorld()
{
	try (var dirStream = Files.walk(Paths.get(GetDirectory()))) {
		dirStream
			.map(Path::toFile)
			.sorted(Comparator.reverseOrder())
			.forEach(File::delete);
	}
	catch (IOException e) {
		e.printStackTrace();
	}
	FileAvailable = false;
}
public void LoadWorld()
{
	LoadPerlinData();
	LoadRegionGridInfo();
	LoadAllSamplePoints();
}*/
void 
WorldMap::EnableAllRegionalRenderings()
{
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i * h + j] != nullptr)
				regions[i * h + j]->EnableRendering();
		}
	}
}
void 
WorldMap::InitNewWorld()
{
	w = 1;
	h = 1;
	regions = new RegionalMap* [w * h];
	x0 = 0;
	y0 = 0;
	RegionalMap::ScrollOriginOffsetForOptimalCoastliness();
	regions[0] = new RegionalMap(x0, y0, this);
	InitPoissonDiscSample(regions[0]);
	//SavePerlinData();
	//SaveRegionGridInfo();
}
/*private boolean LoadAllSamplePoints()
{
	boolean hadProblem = false;
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i][j] == null)
				continue;
			if (!regions[i][j].LoadSampleList())
				hadProblem = true;
		}
	}
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i][j] == null)
				continue;
			if (!regions[i][j].LoadSampleAdjacencies())
				hadProblem = true;
		}
	}
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i][j] == null)
				continue;
			if (!regions[i][j].LoadRiverFile())
				hadProblem = true;
		}
	}
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i][j] == null)
				continue;
			if (!regions[i][j].LoadMidpointTree())
				hadProblem = true;
		}
	}
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (regions[i][j] == null)
				continue;
			regions[i][j].EnableRendering();
		}
	}
	return !hadProblem;
}
private String GetRegionsDescName()
{
	return GetDirectory() + File.separator + "RegionsDesc.txt";
}
private boolean SaveRegionGridInfo()
{
	File regionData = new File(GetRegionsDescName());
	if (regionData.exists())
	{
		if (!regionData.delete())
			return false;
	}
	try {
		BufferedWriter wr = new BufferedWriter(new FileWriter(regionData));
		wr.write(Integer.toString(x0));
		wr.write(" ");
		wr.write(Integer.toString(y0));
		wr.newLine();
		wr.write(Integer.toString(w));
		wr.write(" ");
		wr.write(Integer.toString(h));
		wr.newLine();
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				if (regions[i][j] == null)
					wr.write(K_EMPTY_REGION_NAME);
				else
					wr.write(regions[i][j].GetDirectoryName());
				wr.newLine();
			}
		}
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadRegionGridInfo()
{
	File regionData = new File(GetRegionsDescName());
	if (!regionData.exists() || !regionData.isFile())
		return false;
	try {
		Scanner std = new Scanner(regionData);
		x0 = std.nextInt();
		y0 = std.nextInt();
		std.nextLine();
		w = std.nextInt();
		h = std.nextInt();
		std.nextLine();

		regions = new RegionalMap[w][h];
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				String name = std.nextLine();
				if (name.equals(K_EMPTY_REGION_NAME))
				{
					regions[i][j] = null;
					continue;
				}
				regions[i][j] = new RegionalMap(name, this);
			}
		}
		std.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private String GetPerlinDescName()
{
	return GetDirectory() + File.separator + "PerlinDesc.txt";
}
private boolean SavePerlinData()
{
	File perlinData = new File(GetPerlinDescName());
	if (perlinData.exists())
	{
		if (!perlinData.delete())
			return false;
	}
	try {
		BufferedWriter wr = new BufferedWriter(new FileWriter(perlinData));
		wr.write(Integer.toString(RegionalMap.ORIGIN_OFFSET));
		wr.newLine();
		boolean success = Perlin.SaveSeeds(wr);
		wr.close();
		return success;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadPerlinData()
{
	File perlinData = new File(GetPerlinDescName());
	if (!perlinData.exists() || !perlinData.isFile())
		return false;
	try {
		Scanner std = new Scanner(perlinData);
		RegionalMap.ORIGIN_OFFSET = std.nextInt();
		std.nextLine();
		boolean success = Perlin.LoadSeeds(std);
		std.close();
		return success;
	}
	catch (IOException e) {
		return false;
	}
}*/

RegionalMap* 
WorldMap::GetRegion(int x, int y)
{
	int arrX = x + x0;
	int arrY = y + y0;
	if (arrX >= w || arrY >= h || arrX < 0 || arrY < 0)
		return nullptr;
	return regions[arrX * h + arrY];
}

void 
WorldMap::InitPoissonDiscSample(RegionalMap* target)
{
	vector<SamplePoint*> active;
	SamplePoint* start = new SamplePoint(
		Rand::Double() + target->GetWorldX(),
		Rand::Double() + target->GetWorldY(),
		target);
	if (!target->SetPoint(start))
	{
		delete start;
		return;
	}
	active.push_back(start);
	vector<SamplePoint*> created = PoissonDiscSample(active);
	target->CalculateAllVoronoiAdjacencies();
	//target.SaveSampleList(true);
	//target.SaveSampleAdjacencies(true);
	//TODO we also need to save rivers here, and maybe something with mesh midpoints?
}
vector<SamplePoint*> 
WorldMap::PoissonDiscSample(vector<SamplePoint*>& active)
{
	vector <SamplePoint*> created;
	while (!active.empty())
	{
		int index = (int)(Rand::Double() * active.size());
		SamplePoint* sel = active[index];
		bool placed = false;
		for (int i = 0; i < 30; i++)
		{
			double theta = 2 * numbers::pi * Rand::Double();
			double dm2 = RegionalMap::MIN_VORONOI_DIST * RegionalMap::MIN_VORONOI_DIST;
			double rad = 0;
			if (Switches::POISSON_DENSE)
				rad = (1 + 0.1 * Rand::Double()) * RegionalMap::MIN_VORONOI_DIST;
			else if (Switches::POISSON_BIASED)
				rad = Rand::Double() * RegionalMap::MIN_VORONOI_DIST + RegionalMap::MIN_VORONOI_DIST;
			else
				rad = sqrt(Rand::Double() * (3 * dm2) + dm2);

			double x = sel->x + rad * cos(theta);
			double y = sel->y + rad * sin(theta);
			int indeX = (int)x - RegionalMap::ORIGIN_OFFSET;
			int indeY = (int)y - RegionalMap::ORIGIN_OFFSET;
			RegionalMap* target = GetRegion(indeX, indeY);
			if (target == nullptr)
				continue;
			SamplePoint* cand = new SamplePoint(x, y, target);
			if (target->ExistingPointNearby(cand))
			{
				delete cand;
				cand = nullptr;
				continue;
			}
			if (target->SetPoint(cand))
			{
				created.push_back(cand);
				placed = true;
				active.push_back(cand);
				break;
			}
			else
			{
				delete cand;
				cand = nullptr;
			}
		}
		if (!placed)
		{
			SamplePoint* last = active[active.size() - 1];
			active[index] = last;
			active.pop_back();
		}
	}
	return created;
}
void 
WorldMap::FillAllContinents(int regionX, int regionY, vector<SamplePoint*>& allFound)
{
	deque<SamplePoint*> continentalQueue;
	RegionalMap* target = RequestRegion(regionX, regionY);
	if (target == nullptr)
		return;
	SamplePoint::StartNewSearch();
	for (SamplePoint* p : target->GetAllPoints())
	{
		if (p->type.IsTerrainOfType(TerrainTemplate::OCEAN))
			continue;
		continentalQueue.push_back(p);
		p->MarkAsReached();
	}
	while (!continentalQueue.empty())
	{
		SamplePoint* v = continentalQueue.front();
		continentalQueue.pop_front();
		allFound.push_back(v);
		if (v->NearNorthEdge() && v->GetRegionalMap()->GetNorth() == nullptr)
		{
			RequestRegion(v->GetRegionalMap()->GetOriginX(), v->GetRegionalMap()->GetOriginY() - 1);
		}
		if (v->NearSouthEdge() && v->GetRegionalMap()->GetSouth() == nullptr)
		{
			RequestRegion(v->GetRegionalMap()->GetOriginX(), v->GetRegionalMap()->GetOriginY() + 1);
		}
		if (v->NearWestEdge() && v->GetRegionalMap()->GetWest() == nullptr)
		{
			RequestRegion(v->GetRegionalMap()->GetOriginX() - 1, v->GetRegionalMap()->GetOriginY());
		}
		if (v->NearEastEdge() && v->GetRegionalMap()->GetEast() == nullptr)
		{
			RequestRegion(v->GetRegionalMap()->GetOriginX() + 1, v->GetRegionalMap()->GetOriginY());
		}
		for (SamplePoint* adj : v->GetAdjacentSamples())
		{
			if (adj->Reached())
				continue;
			if (adj->type.IsTerrainOfType(TerrainTemplate::OCEAN))
				continue;
			adj->MarkAsReached();
			continentalQueue.push_back(adj);
		}
	}
}
RegionalMap* 
WorldMap::RequestRegion(int x, int y)
{
	int arrX = x + x0;
	int arrY = y + y0;
	if (arrX < w && arrY < h && arrX >= 0 && arrY >= 0)
	{
		if (regions[arrX * h + arrY] != nullptr)
			return regions[arrX * h + arrY];
		CreateRegion(x, y);
		return GetRegion(x, y);
	}
	int posXInc = arrX - w + 1;
	if (posXInc < 0)
		posXInc = 0;
	int negXInc = arrX * -1;
	if (negXInc < 0)
		negXInc = 0;
	int posYInc = arrY - h + 1;
	if (posYInc < 0)
		posYInc = 0;
	int negYInc = arrY * -1;
	if (negYInc < 0)
		negYInc = 0;
	ExpandEmptySpace(posXInc, negXInc, posYInc, negYInc);
	CreateRegion(x, y);
	return GetRegion(x, y);
}
bool 
WorldMap::CreateRegion(int x, int y)
{
	int arrX = x + x0;
	int arrY = y + y0;
	if (arrX >= w || arrY >= h || arrX < 0 || arrY < 0)
		return false;
	if (regions[arrX * h + arrY] != nullptr)
		return false;

	regions[arrX * h + arrY] = new RegionalMap(x, y, this);
	vector<SamplePoint*> newlyActive;
	RegionalMap* north = GetRegion(x, y - 1);
	if (north != nullptr)
	{
		for (int ii = 0; ii < RegionalMap::VORONOI_DIM; ii++)
		{
			SamplePoint* p1 = north->GetAt(ii, RegionalMap::VORONOI_DIM - 1);
			SamplePoint* p2 = north->GetAt(ii, RegionalMap::VORONOI_DIM - 2);
			if (p1 != nullptr)
				newlyActive.push_back(p1);
			if (p2 != nullptr)
				newlyActive.push_back(p2);
		}
	}
	RegionalMap* south = GetRegion(x, y + 1);
	if (south != nullptr)
	{
		for (int ii = 0; ii < RegionalMap::VORONOI_DIM; ii++)
		{
			SamplePoint* p1 = south->GetAt(ii, 0);
			SamplePoint* p2 = south->GetAt(ii, 1);
			if (p1 != nullptr)
				newlyActive.push_back(p1);
			if (p2 != nullptr)
				newlyActive.push_back(p2);
		}
	}
	RegionalMap* west = GetRegion(x - 1, y);
	if (west != nullptr)
	{
		for (int jj = 0; jj < RegionalMap::VORONOI_DIM; jj++)
		{
			SamplePoint* p1 = west->GetAt(RegionalMap::VORONOI_DIM - 1, jj);
			SamplePoint* p2 = west->GetAt(RegionalMap::VORONOI_DIM - 2, jj);
			if (p1 != nullptr)
				newlyActive.push_back(p1);
			if (p2 != nullptr)
				newlyActive.push_back(p2);
		}
	}
	RegionalMap* east = GetRegion(x + 1, y);
	if (east != nullptr)
	{
		for (int jj = 0; jj < RegionalMap::VORONOI_DIM; jj++)
		{
			SamplePoint* p1 = east->GetAt(0, jj);
			SamplePoint* p2 = east->GetAt(1, jj);
			if (p1 != nullptr)
				newlyActive.push_back(p1);
			if (p2 != nullptr)
				newlyActive.push_back(p2);
		}
	}
	if (newlyActive.size() == 0)
	{
		InitPoissonDiscSample(regions[arrX * h + arrY]);
	}
	else
	{
		vector<SamplePoint*> newV = PoissonDiscSample(newlyActive);
		unordered_set<SamplePoint*> newlyCreated;
		unordered_set<SamplePoint*> corrupted;
		for (SamplePoint* v : newV)
		{
			v->CalculateAdjacencies();
		}
		for (SamplePoint* v : newV)
		{
			for (SamplePoint* vv : v->GetAdjacentSamples())
			{
				if (!newlyCreated.contains(vv))
					corrupted.insert(vv);
			}
		}
		for (SamplePoint* v : corrupted)
			v->ResetAdjacencies();
		for (SamplePoint* v : corrupted)
			v->CalculateAdjacencies();

		//SaveSampleFiles(newV, true, true, false, false);
		//SaveSampleFiles(corrupted, false, true, false, false);
	}
	//SaveRegionGridInfo();
	return true;
}

void 
WorldMap::ExpandEmptySpace(int posXInc, int negXInc, int posYInc, int negYInc)
{
	int newW = w + posXInc + negXInc;
	int newH = h + posYInc + negYInc;
	int newX0 = x0 + negXInc;
	int newY0 = y0 + negYInc;
	RegionalMap** newRegions = new RegionalMap * [newW * newH] {};
	for (int i = 0; i < newW; i++)
		for (int j = 0; j < newH; j++)
		{
			int oldX = i - negXInc;
			int oldY = j - negYInc;
			if (oldX >= 0 && oldY >= 0 && oldX < w && oldY < h)
			{
				newRegions[i * newH + j] = regions[oldX * h + oldY];
			}
		}
	w = newW;
	h = newH;
	x0 = newX0;
	y0 = newY0;
	delete regions;
	regions = newRegions;
	//SaveRegionGridInfo();
}
/*public void SaveSampleFiles(Collection<SamplePoint> points,
	boolean basicData,
	boolean adjacencies,
	boolean midpointTrees,
	boolean rivers)
{
	HashSet<RegionalMap> updated = new HashSet<RegionalMap>();
	for (SamplePoint v : points)
	{
		updated.add(v.GetRegionalMap());
	}
	for (RegionalMap rm : updated)
	{
		if (basicData)
			rm.SaveSampleList(true);
		if (adjacencies)
			rm.SaveSampleAdjacencies(true);
		if (midpointTrees)
			rm.SaveMidpointTree(true);
		if (rivers)
			rm.SaveRiverFile(true);
	}
}
public void AttachTool(WorldMapTool newTool)
{
	DettachActiveTool();
	if (newTool != null)
	{
		if (!newTool.ZoomOkay())
			StopListeningToZoom();
		if (!newTool.DragOkay())
			StopListeningToDrag();
		newTool.Activate();
		activeTool = newTool;
		addMouseListener(activeTool);
		addMouseMotionListener(activeTool);
		addKeyListener(activeTool);
		this.requestFocusInWindow();
	}
}
public void DettachActiveTool()
{
	if (activeTool != null)
	{
		activeTool.Deactivate();
		removeMouseListener(activeTool);
		removeMouseMotionListener(activeTool);
		removeKeyListener(activeTool);
		if (!activeTool.DragOkay())
			StartListeningToDrag();
		if (!activeTool.ZoomOkay())
			StartListeningToZoom();
	}
}*/
RegionalMap::Coordinate 
WorldMap::GetRegionalMapAt(double x, double y)
{
	return GetRegionalMapAt(x, y, false);
}
RegionalMap::Coordinate 
WorldMap::GetRegionalMapAt(double x, double y, bool createIfNot)
{
	double regionDim = tileSize * RegionalMap::DIMENSION;
	//x -= getWidth() / 2;
	//y -= getHeight() / 2;
	x += x0 * regionDim;
	y += y0 * regionDim;
	x += 0.5 * regionDim;
	y += 0.5 * regionDim;
	x -= dX;
	y -= dY;

	int xIndex = (int) floor(x / regionDim);
	int yIndex = (int) floor(y / regionDim);

	RegionalMap* found = nullptr;
	if (createIfNot)
		found = RequestRegion(xIndex - x0, yIndex - y0);
	else
		found = GetRegion(xIndex - x0, yIndex - y0);

	if (found == nullptr)
		return RegionalMap::Coordinate(0, 0, nullptr);

	x -= xIndex * regionDim;
	y -= yIndex * regionDim;
	x /= regionDim;
	y /= regionDim;

	return RegionalMap::Coordinate(x, y, found);
}

void
WorldMap::RenderVoronoiTiles(int width, int height, float regionDim, bool terrainBased)
{
	float left = 0.0f;
	float right = (float)width;
	float bottom = 0.0f;
	float top = (float)height;
	float nearZ = -1.0f;
	float farZ = 1.0f;

	float mvdScreenSize = (float) RegionalMap::MIN_VORONOI_DIST * regionDim;

	glm::mat4 orthoProj = glm::ortho(left, right, bottom, top, nearZ, farZ);
	if (!voronoiRenderer.IsInitialized())
		voronoiRenderer.Init();

	vector<glm::vec2> samplePointLocs;
	vector<glm::vec3> samplePointColors;
	ForEachOnSceenRegionalMap(
		[&samplePointLocs, &samplePointColors, terrainBased, regionDim, width, height]
		(RegionalMap* region, float xLoc, float yLoc)
		{
			region->GatherSamplePointCenters(regionDim, width, height, xLoc, yLoc, samplePointLocs);
			region->GatherSamplePointColors(regionDim, width, height, xLoc, yLoc, samplePointColors, terrainBased);
		}, width, height);
	voronoiRenderer.Render(samplePointLocs, samplePointColors, mvdScreenSize * 3, orthoProj);

}


void
WorldMap::RenderRegionalData(int width, int height, float regionDim, RegionalDataRenderMode mode)
{
	float left = 0.0f;
	float right = (float)width;
	float bottom = 0.0f;
	float top = (float)height;
	float nearZ = -1.0f;
	float farZ = 1.0f;

	glm::mat4 orthoProj = glm::ortho(left, right, bottom, top, nearZ, farZ);

	if (!regionalDataRenderer.IsInitialized())
		regionalDataRenderer.Init();

	myWorkerThread.DiscardAllWaitingTasks();

	vector<glm::vec2> sectionLocs;
	vector<RegionalDataLoc> sections;
	vector<RegionalDataCheck*> checks;
	checks.push_back(&heightmaps);

	ForEachOnSceenRegionalMap(
		[this, orthoProj, regionDim, width, height, &sectionLocs, &sections, &checks](RegionalMap* region, float xLoc, float yLoc)
		{
			region->GatherRelevantDataIDs(regionDim, width, height, xLoc, yLoc, checks, sectionLocs, sections, true);
		}, width, height);

	for (int i = 0; i < sections.size(); i++)
	{
		Heightmap hm = heightmaps.GetRaster(sections[i]);
		glm::vec2 renderLoc = sectionLocs[i];
		float dim = regionDim / sections[i].GetNumSectionsPerSide();
		float meterDim = RegionalMap::METER_DIM / sections[i].GetNumSectionsPerSide();
		float mpp = meterDim / hm.DIM;

		glm::vec4 color = glm::vec4(1, 1, 1, 1);
		if (sections[i].LOD == 2)
			color = glm::vec4(1, 0, 0, 1);
		else if (sections[i].LOD == 3)
			color = glm::vec4(0, 1, 0, 1);
		else if(sections[i].LOD == 4)
			color = glm::vec4(0, 0, 1, 1);
		else if (sections[i].LOD == 5)
			color = glm::vec4(1, 1, 0, 1);
		if (hm.OkayToUse())
		{
			switch (mode)
			{
			case BANDED_HEIGHTS:
				regionalDataRenderer.RenderBandedHeight(renderLoc, dim, dim, orthoProj, hm.GetRawData());
				break;
			case HILLSHADE:
				regionalDataRenderer.RenderHillshade(renderLoc, dim, dim, orthoProj, hm.GetRawData(), mpp);
				break;
			}

		}
			
	}

	if (mode != HILLSHADE)
		return;

	sectionLocs.clear();
	sections.clear();
	checks.clear();
	checks.push_back(&depthmaps);

	ForEachOnSceenRegionalMap(
		[this, orthoProj, regionDim, width, height, &sectionLocs, &sections, &checks](RegionalMap* region, float xLoc, float yLoc)
		{
			region->GatherRelevantDataIDs(regionDim, width, height, xLoc, yLoc, checks, sectionLocs, sections, false);
		}, width, height);

	for (int i = 0; i < sections.size(); i++)
	{
		Depthmap dm = depthmaps.GetRaster(sections[i]);
		glm::vec2 renderLoc = sectionLocs[i];
		float dim = regionDim / sections[i].GetNumSectionsPerSide();
		float meterDim = RegionalMap::METER_DIM / sections[i].GetNumSectionsPerSide();
		float mpp = meterDim / dm.DIM;

		glm::vec4 color = glm::vec4(1, 1, 1, 1);
		if (sections[i].LOD == 2)
			color = glm::vec4(1, 0, 0, 1);
		else if (sections[i].LOD == 3)
			color = glm::vec4(0, 1, 0, 1);
		else if (sections[i].LOD == 4)
			color = glm::vec4(0, 0, 1, 1);
		else if (sections[i].LOD == 5)
			color = glm::vec4(1, 1, 0, 1);
		if (!dm.OkayToUse())
			continue;
		regionalDataRenderer.RenderRivers(renderLoc, dim, dim, orthoProj, dm.GetRawData(), mpp);
	}
}

void
WorldMap::RenderMeshPoints(int width, int height, float regionDim)
{
	float left = 0.0f;
	float right = (float)width;
	float bottom = 0.0f;
	float top = (float)height;
	float nearZ = -1.0f;
	float farZ = 1.0f;

	double mvdScreenSize = RegionalMap::MIN_VORONOI_DIST * regionDim;

	if (mvdScreenSize < 12)
		return;

	glm::mat4 orthoProj = glm::ortho(left, right, bottom, top, nearZ, farZ);
	if (!meshPointRenderer.IsInitialized())
		meshPointRenderer.Init();

	vector<glm::vec2> samplePointLocs;
	ForEachOnSceenRegionalMap(
		[&samplePointLocs, regionDim, width, height](RegionalMap* region, float xLoc, float yLoc)
		{
			region->GatherSamplePointCenters(regionDim, width, height, xLoc, yLoc, samplePointLocs);
		}, width, height);
	meshPointRenderer.Render(samplePointLocs, (float) 6, (float) 0, orthoProj, glm::vec4(0, 0, 0, 1));

	if (mvdScreenSize < 24)
		return;
	meshPointRenderer.Render(samplePointLocs, (float)mvdScreenSize, (float)mvdScreenSize - 2, orthoProj, 
		glm::vec4(1.0f, 186.f/255.f, 0.0f, 1));

}

void 
WorldMap::RenderOutlines(int width, int height, float regionDim)
{

	float left = 0.0f;
	float right = (float) width;
	float bottom = 0.0f;
	float top = (float) height;
	float nearZ = -1.0f;
	float farZ = 1.0f;

	glm::mat4 orthoProj = glm::ortho(left, right, bottom, top, nearZ, farZ);
	if (!outlineRenderer.IsInitialized())
		outlineRenderer.Init();

	int lmSize = this->tileSize;
	vector<glm::vec2> outlineLocs;
	ForEachOnSceenRegionalMap(
		[&outlineLocs, regionDim, lmSize, width, height](RegionalMap* region, float xLoc, float yLoc)
		{
			region->GatherLocalMapOutlineLocations(lmSize, width, height, xLoc, yLoc, outlineLocs);
		}, width, height);
	outlineRenderer.Render(outlineLocs, (float) tileSize, (float) tileSize, orthoProj, glm::vec4(1, 0, 0, 1));

	outlineLocs.clear();
	ForEachOnSceenRegionalMap(
		[&outlineLocs, regionDim, lmSize, width, height](RegionalMap* region, float xLoc, float yLoc)
		{
			outlineLocs.push_back(glm::vec2(xLoc, yLoc));
		}, width, height);
	outlineRenderer.Render(outlineLocs, regionDim, regionDim, orthoProj, glm::vec4(0, 1, 0, 1));
}

GLuint
WorldMap::Render(int width, int height) 
{
	if(!textureInit || texWidth != width || texHeight != height) {
		InitializeRenderTarget(width, height);
		textureInit = true;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
	glViewport(0, 0, width, height);  // Set to texture size

	glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	float regionDim = (float)(tileSize * RegionalMap::DIMENSION);
	// Render your map here...
	switch (Switches::CURR_PAINT_TYPE)
	{
	case Switches::PAINT_TYPE::TERRAIN:
		RenderVoronoiTiles(width, height, regionDim, true);
		break;
	case Switches::PAINT_TYPE::VORONOI_PURE:
		RenderVoronoiTiles(width, height, regionDim, false);
		break;
	case Switches::PAINT_TYPE::ELEVATION_CURR:
		RenderRegionalData(width, height, regionDim, RegionalDataRenderMode::BANDED_HEIGHTS);
		break;
	case Switches::PAINT_TYPE::ELEV_GRADIENT:
		RenderRegionalData(width, height, regionDim, RegionalDataRenderMode::HILLSHADE);
		break;
	case Switches::PAINT_TYPE::CONTOUR:
		RenderRegionalData(width, height, regionDim, RegionalDataRenderMode::CONTOUR);
		break;
	default:
		break;
	}
	if (Switches::OUTLINE_MAPS)
	{
		RenderOutlines(width, height, regionDim);
	}
	if (Switches::PAINT_VORONOI_CENTERS)
	{
		RenderMeshPoints(width, height, regionDim);
	}
	//myWorkerThread.discardAllQueuedTasks();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);  // Return to default framebuffer

	return renderTexture;  // Return the texture ID for use in rendering
}

void 
WorldMap::InitializeRenderTarget(int width, int height) 
{

	if (frameBuffer != 0) {
		glDeleteFramebuffers(1, &frameBuffer);
		glDeleteTextures(1, &renderTexture);
	}

	texWidth = width;
	texHeight = height;

	glGenFramebuffers(1, &frameBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

	// Create texture to render to
	glGenTextures(1, &renderTexture);
	glBindTexture(GL_TEXTURE_2D, renderTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D, renderTexture, 0);

	glGenTextures(1, &depthBuffer);
	glBindTexture(GL_TEXTURE_2D, depthBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthBuffer, 0);
}
