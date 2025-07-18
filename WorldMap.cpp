#include "WorldMap.h"
#include "RegionalMap.h"
#include "SamplePoint.h"
#include "GlobalRand.h"
#include "Switches.h"
#include <cmath>
#include <numbers>
#include <deque>
#include <unordered_set>

using namespace std;
/*public WorldMap(String name)
{
	this.worldName = name;
	tileSize = DEFAULT_TILE_SIZE;
	File topDir = new File(GetDirectory());
	if (topDir.exists() && topDir.isDirectory())
	{
		FileAvailable = true;
	}
	else
	{
		FileAvailable = false;
		topDir.mkdir();
		InitNewWorld();
	}
}
public String GetDirectory()
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
	regions = new RegionalMap * [w * h];
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
	PoissonDiscSample(active);
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
	RegionalMap** newRegions = new RegionalMap * [newW * newH];
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
}
public void StartListeningToZoom()
{
	StopListeningToZoom();
	myZoom = new ZoomUpdater();
	addMouseWheelListener(myZoom);
}
public void StopListeningToZoom()
{
	if (myZoom != null)
	{
		removeMouseWheelListener(myZoom);
		myZoom = null;
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
	//x -= dX;
	//y -= dY;

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
public synchronized void ZoomIn(int mouseX, int mouseY)
{
	double proposedZoom = 1.1 * tileSize;
	if (proposedZoom > MAXIMUM_TILE_SIZE)
		return;
	int actualZoom = (int)(proposedZoom + 0.5);
	double actualRatio = 1.0 * actualZoom / tileSize;
	tileSize = actualZoom;
	AdjustDeltas(mouseX, mouseY, actualRatio);
}
public synchronized void ZoomOut(int mouseX, int mouseY)
{
	double proposedZoom = tileSize / 1.1;
	if (proposedZoom < MINIMUM_TILE_SIZE)
		return;
	int actualZoom = (int)proposedZoom;
	double actualRatio = 1.0 * actualZoom / tileSize;
	tileSize = actualZoom;
	AdjustDeltas(mouseX, mouseY, actualRatio);
}
public double GetTileSize()
{
	return tileSize;
}
private void AdjustDeltas(int mouseX, int mouseY, double ratio)
{
	mouseX -= getWidth() / 2;
	mouseY -= getHeight() / 2;
	double newDX = dX * ratio - mouseX * (ratio - 1);
	double newDY = dY * ratio - mouseY * (ratio - 1);
	dX = newDX;
	dY = newDY;
}
public synchronized long Render()
{
	BufferedImage buffer = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
	Graphics2D g2 = buffer.createGraphics();
	double regionDim = tileSize * RegionalMap.DIMENSION;
	g2.setColor(Color.black);
	g2.fillRect(0, 0, getWidth(), getHeight());
	g2.setClip(0, 0, getWidth(), getHeight());
	AffineTransform saved = g2.getTransform();
	g2.translate(getWidth() / 2, getHeight() / 2);
	g2.translate(-1 * x0 * regionDim, -1 * y0 * regionDim);
	g2.translate(-0.5 * regionDim, -0.5 * regionDim);
	g2.translate(dX, dY);
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			if (g2.hitClip(0, 0, (int)regionDim, (int)regionDim))
			{
				RegionalMap m = regions[i][j];
				if (m != null)
					m.Render(regionDim, g2);
			}

			g2.translate(0, regionDim);
		}
		g2.translate(regionDim, -1 * h * regionDim);
	}
	RenderQueue.IncrementFrameID();
	g2.setTransform(saved);
	RenderScale(g2);
	Graphics g = getGraphics();
	if (g == null)
		return 10l;
	g.drawImage(buffer, 0, 0, null);
	return 10l;
}
private void RenderScale(Graphics2D g2)
{
	double mPerTile = LocalMap.METER_DIM;
	double mPerPixel = mPerTile / tileSize;
	double mInTargetScale = SCALE_LABEL_TARGET_WIDTH * mPerPixel;
	double mInChosenScale = -1;
	String chosenScale = "";
	for (int i = 0; i < LABEL_DISTANCES.length; i++)
	{
		double mInScale = LABEL_DISTANCES[i];
		String scaleName = TARGET_SCALE_LABLES[i];
		if (i == 0 || (Math.abs(mInTargetScale - mInScale) < Math.abs(mInTargetScale - mInChosenScale)))
		{
			mInChosenScale = mInScale;
			chosenScale = scaleName;
		}
	}
	Font f = new Font("Calibri", Font.BOLD, 18);
	//Font f = g2.getFont().deriveFont(Font.BOLD, 16);
	g2.setFont(f);
	double barWidth = mInChosenScale / mPerPixel;
	int textWidth = g2.getFontMetrics().stringWidth(chosenScale);

	int h = g2.getFontMetrics().getHeight() + 15;
	int w = (int)Math.max(textWidth, barWidth) + 10;

	AffineTransform saved = g2.getTransform();
	g2.translate(getWidth(), getHeight());
	g2.translate(-1 * w, -1 * h);
	g2.translate(-10, -10);
	g2.translate(w / 2, 0);
	g2.setColor(Color.white);
	g2.fillRect((int)(-0.5 * textWidth - 5), 0, textWidth + 10, g2.getFontMetrics().getHeight());
	g2.setColor(Color.black);
	g2.drawString(chosenScale, (int)(-0.5 * textWidth), (int)(0.75 * g2.getFontMetrics().getHeight()));
	g2.translate(0, h / 2);
	g2.setColor(Color.orange);
	int lineY = (int)(h / 2 - 5);
	g2.fillRect((int)(-0.5 * barWidth), lineY - 2, (int)(barWidth), 4);
	g2.fillRect((int)(-0.5 * barWidth - 1), lineY - 6, 2, 12);
	g2.fillRect(-1, lineY - 4, 2, 8);
	g2.fillRect((int)(0.5 * barWidth - 1), lineY - 6, 2, 12);
	g2.setTransform(saved);
}
private class ZoomUpdater implements MouseWheelListener
{

	@Override
		public void mouseWheelMoved(MouseWheelEvent e) {
		int rots = e.getWheelRotation();
		if (rots < 0)
			for (int i = 0; i < rots * -1; i++)
				ZoomIn(e.getX(), e.getY());
		else
			for (int i = 0; i < rots; i++)
				ZoomOut(e.getX(), e.getY());
	}
}
