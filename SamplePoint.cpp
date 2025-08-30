#include "SamplePoint.h"
#include "GlobalRand.h"
#include "Perlin.h"
#include "MeshConnection.h"
#include "RegionalMap.h"

using namespace std;


SamplePoint::SamplePoint(double x, double y, RegionalMap* parent) : 
	MeshPoint(x, y), paintMyColor(false), parent(parent), riverFlow(0), 
	riverOutlet(nullptr), maxGrade(0), tectonicUpliftOverride(0)
{
	float r = Rand::Float();
	float g = Rand::Float();
	float b = Rand::Float();
	myColor = glm::vec3(r, g, b);
	Init();
}

void 
SamplePoint::Init()
{
	paintMyColor = false;
	maxGrade = 0;
	ResetRiver();
	ResetTectonicUplift();
}
/*public SamplePoint(DataInputStream dis, RegionalMap parent)
{
	super(dis);
	try
	{
		short typeBits = dis.readShort();
		int rgb = dis.readInt();
		type = new TerrainType(typeBits);
		myColor = new Color(rgb);
		this.parent = parent;
		Init();
	}
	catch (IOException e)
	{
		e.printStackTrace();
	}
}
public SamplePoint(Iterator<String> tokenStream, RegionalMap parent)
{
	super(tokenStream);
	short typeBits = Short.parseShort(tokenStream.next());
	int rgb = Integer.parseInt(tokenStream.next());
	type = new TerrainType(typeBits);
	myColor = new Color(rgb);
	this.parent = parent;
	Init();
}
public void SetAdjacenciesFromDescription(String desc, WorldMap wm)
{
	String[] tokens = desc.split(" ");
	for (int i = 0; i < tokens.length; i++)
	{
		SamplePoint adj = FindFromGlobalIdentifier(tokens[i], wm);
		if (adj == null)
			continue;
		MarkAdjacent(adj);
	}
}
public static SamplePoint FindFromGlobalIdentifier(String desc, WorldMap wm)
{
	String[] dcr = desc.split(",");
	int px = Integer.parseInt(dcr[0]);
	int py = Integer.parseInt(dcr[1]);
	int index = Integer.parseInt(dcr[2]);
	RegionalMap par = wm.GetRegion(px, py);
	if (par == null)
		return null;
	return par.GetSamplePoint(index);
}
public boolean LoadAdjacencies(DataInputStream dis, WorldMap wm)
{
	try {
		int num = dis.readByte();
		for (int i = 0; i < num; i++)
		{
			SamplePoint adj = LoadFromGlobalIdentifier(dis, wm);
			MarkAdjacent(adj);
		}
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public static SamplePoint LoadFromGlobalIdentifier(DataInputStream dis, WorldMap wm)
{
	try
	{
		int px = dis.readByte();
		int py = dis.readByte();
		int index = dis.readInt();
		RegionalMap par = wm.GetRegion(px, py);
		if (par == null)
		{
			return null;

		}
		return par.GetSamplePoint(index);
	}
	catch (IOException e)
	{
		return null;
	}
}
public boolean WriteAdjacencies(DataOutputStream dos)
{
	try
	{
		boolean worked = true;
		dos.writeByte(GetAdjacentSamples().size());
		for (SamplePoint adj : GetAdjacentSamples())
		{
			worked = worked && adj.WriteGlobalIdentifier(dos);
		}
		return worked;
	}
	catch (IOException e)
	{
		return false;
	}
}
public String GetAdjacencyDescription()
{
	String desc = "";
	ArrayList<SamplePoint> adjacent = GetAdjacentSamples();
	for (int i = 0; i < adjacent.size(); i++)
	{
		SamplePoint adj = adjacent.get(i);
		desc += adj.GetGlobalIdentifier();
		if (i != adjacent.size() - 1)
			desc += " ";
	}
	return desc;
}
public String GetGlobalIdentifier()
{
	String desc = Integer.toString(parent.GetOriginX()) + ",";
	desc += Integer.toString(parent.GetOriginY()) + ",";
	desc += Integer.toString(GetContainerIndex());
	return desc;
}
public boolean WriteGlobalIdentifier(DataOutputStream dos)
{
	try
	{
		dos.writeByte(parent.GetOriginX());
		dos.writeByte(parent.GetOriginY());
		dos.writeInt(GetContainerIndex());
		return true;
	}
	catch (IOException e)
	{
		return false;
	}
}
public boolean WriteDescription(DataOutputStream dos)
{
	if (!super.WriteDescription(dos))
		return false;
	try {
		dos.writeShort(type.BitsToShort());
		dos.writeInt(myColor.getRGB());
		return true;
	}
	catch (IOException e) {
		return false;
	}
}

public String GetDescription()
{
	String desc = super.GetDescription() + " ";
	desc += Short.toString(type.BitsToShort()) + " ";
	desc += Integer.toString(myColor.getRGB());
	return desc;
}*/

bool 
SamplePoint::RidgeAssigned()
{
	bool result = false;
	ForEachAdjacent([this, &result](MeshPoint* p) {
		MeshConnection* m = GetConnection(p);
		if (m->IsRidgeline())
			result = true;
	});
	return result;
}

bool
SamplePoint::SetRiverOutlet(SamplePoint* p)
{
	MeshConnection* m = GetConnection(p);
	if (m == nullptr)
		return false;
	if (riverOutlet != nullptr)
	{
		MeshConnection* res = GetConnection(riverOutlet);
		if (res == nullptr)
			return false;
		res->ResetRiver();
	}
	m->SetRiver();
	riverOutlet = p;
	p->riverInlets.push_back(this);
	return true;
}
void
SamplePoint::SendFlowToOutlet()
{
	if (riverOutlet != nullptr)
		riverOutlet->riverFlow += riverFlow;
}
vector<SamplePoint*>
SamplePoint::GrabAdjacentUnassignedRiverFlows(bool sendRiverFlow)
{
	vector<SamplePoint*> grabbed;
	for (SamplePoint* a : GetAdjacentSamples())
	{
		if (a->GetElevation() != GetElevation())
			continue;
		if (a->riverOutlet != nullptr)
			continue;

		a->SetRiverOutlet(this);
		grabbed.push_back(a);
		if (sendRiverFlow)
			a->SendFlowToOutlet();
	}
	return grabbed;
}
SamplePoint* 
SamplePoint::GetWayDownhill(bool markOutlet, bool sendRiverFlow)
{
	double maxDelta = -1;
	SamplePoint* lowest = nullptr;
	for (SamplePoint* p : GetAdjacentSamples())
	{
		double delta = GetElevation() - p->GetElevation();
		if (delta > 0 && delta > maxDelta)
		{
			lowest = p;
			maxDelta = delta;
		}
	}
	if (markOutlet && lowest != nullptr)
	{
		SetRiverOutlet(lowest);
		if (sendRiverFlow)
			lowest->riverFlow += riverFlow;
	}
	return lowest;
}
void 
SamplePoint::ResetRiver()
{
	riverFlow = 1;
	if (riverOutlet != nullptr)
	{
		MeshConnection* m = GetConnection(riverOutlet);
		if (m != nullptr)
			m->ResetRiver();
	}
	riverOutlet = nullptr;
	riverInlets.clear();
}

void 
SamplePoint::SendRiverFlow()
{
	SamplePoint* lowest = GetWayDownhill(true, true);
	if (lowest == nullptr)
	{
		type.ApplyTerrain(TerrainTemplate::LAKE);
	}
}

bool 
SamplePoint::NearNorthEdge()
{
	int j = (int)((y - parent->GetWorldY()) * RegionalMap::VORONOI_DIM);
	return j < 3;
}
bool 
SamplePoint::NearSouthEdge()
{
	int j = (int)((y - parent->GetWorldY()) * RegionalMap::VORONOI_DIM);
	return j >= RegionalMap::VORONOI_DIM - 3;
}
bool 
SamplePoint::NearWestEdge()
{
	int i = (int)((x - parent->GetWorldX()) * RegionalMap::VORONOI_DIM);
	return i < 3;
}
bool 
SamplePoint::NearEastEdge()
{
	int i = (int)((x - parent->GetWorldX()) * RegionalMap::VORONOI_DIM);
	return i >= RegionalMap::VORONOI_DIM - 3;
}

void 
SamplePoint::CalculateAdjacencies()
{
	//10 seems sufficient to find all adjacencies, but it occasionally messes up
	//the circumcenter check on real adjacencies. Boosting to 15 seems sufficient
	//to make that check effectively perfect
	const vector<SamplePoint*> adj = parent->GetNearestN(x, y, 15);
	for (SamplePoint* vp : adj)
	{
		double dist = sqrt((vp->x - x) * (vp->x - x) + (vp->y - y) * (vp->y - y));
		dist /= 2;

		double mx = (vp->x + x) / 2;
		double my = (vp->y + y) / 2;

		vector<SamplePoint*> objectors;
		for (int i = 0; i < adj.size(); i++)
		{
			SamplePoint* a = adj[i];
			if (a == vp)
				continue;
			double compDist = sqrt((a->x - mx) * (a->x - mx) + (a->y - my) * (a->y - my));
			if (compDist < dist)
				objectors.push_back(a);
		}
		if (objectors.empty())
		{
			MarkAdjacent(vp);
			continue;
		}
		bool objection = false;
		for (SamplePoint* b : objectors)
		{
			//now we need to find the circumcenter of this, vp, and b; I'm taking this calc from Wikipedia
			double D = 2 * (x * (vp->y - b->y) + vp->x * (b->y - y) + b->x * (y - vp->y));
			double ccx = ((x * x + y * y) * (vp->y - b->y) + 
				(vp->x * vp->x + vp->y * vp->y) * (b->y - y) + 
				(b->x * b->x + b->y * b->y) * (y - vp->y)) / D;
			double ccy = ((x * x + y * y) * (b->x - vp->x) + 
				(vp->x * vp->x + vp->y * vp->y) * (x - b->x) + 
				(b->x * b->x + b->y * b->y) * (vp->x - x)) / D;

			double ccr = sqrt((x - ccx) * (x - ccx) + (y - ccy) * (y - ccy));
			//double ccr2 = Math.sqrt((vp.x - ccx) * (vp.x - ccx) + (vp.y - ccy) * (vp.y - ccy));
			//double ccr3 = Math.sqrt((b.x - ccx) * (b.x - ccx) + (b.y - ccy) * (b.y - ccy));

			for (int i = 0; i < adj.size(); i++)
			{
				SamplePoint* a = adj[i];
				if (a == vp)
					continue;
				if (a == b)
					continue;
				double compDist = sqrt((a->x - ccx) * (a->x - ccx) + (a->y - ccy) * (a->y - ccy));
				if (compDist < ccr)
				{
					objection = true;
					break;
				}
			}
		}

		if (!objection)
		{
			MarkAdjacent(vp);
		}
	}
}

vector<SamplePoint*> 
SamplePoint::GetAdjacentSamples()
{
	vector<SamplePoint*> samples;
	ForEachAdjacent([&samples](MeshPoint* a) {
		SamplePoint* s = dynamic_cast<SamplePoint*>(a);
		if (s != nullptr)
			samples.push_back(s);
	});
	return samples;
}
void 
SamplePoint::MakeLake()
{
	if (type.IsTerrainOfType(TerrainTemplate::OCEAN))
		type = TerrainType(TerrainTemplate::FLAT);
	type.ApplyTerrain(TerrainTemplate::LAKE);
}
void 
SamplePoint::MakeOcean()
{
	SetElevation(0);
	type = TerrainType(TerrainTemplate::OCEAN);
}

void 
SamplePoint::AssignMaxGrade(double min, double max)
{
	double perlinPush = Perlin::minMaxSelector.Get(x, y);
	perlinPush += 1;
	perlinPush /= 2;
	if (perlinPush < 0)
		perlinPush = 0;
	if (perlinPush > 1)
		perlinPush = 1;
	double grade = perlinPush * max + (1 - perlinPush) * min;
	maxGrade = grade;
}
void 
SamplePoint::SetPeaks()
{
	type = TerrainType(TerrainTemplate::PEAKS);
	AssignMaxGrade(.5, .940);
}
void 
SamplePoint::SetMountains()
{
	type = TerrainType(TerrainTemplate::MOUNTAINS);
	AssignMaxGrade(.342, .866);
}
void 
SamplePoint::SetOcean()
{
	type = TerrainType(TerrainTemplate::OCEAN);
	SetElevation(0);
	AssignMaxGrade(0, 0);
}
void 
SamplePoint::SetHills()
{
	type = TerrainType(TerrainTemplate::HILLS);
	AssignMaxGrade(.174, .5);
}
void 
SamplePoint::SetLakes()
{
	type = TerrainType(TerrainTemplate::FLAT);
	type.ApplyTerrain(TerrainTemplate::LAKE);
	AssignMaxGrade(0, 0);
}
void 
SamplePoint::SetFlats()
{
	TerrainType def = TerrainType(TerrainTemplate::FLAT);
	AssignMaxGrade(.017, .139);
	type = def;
}
void 
SamplePoint::AssignVoronoiTerrainType()
{
	if (!Perlin::oceans.UnderThreshold(x, y))
	{
		if (Perlin::peaks.UnderThreshold(x, y) &&
			Perlin::peaks.GetPercentBeneathThreshold(x, y) > 0.1 &&
			Perlin::oceans.GetPercentAboveThreshold(x, y) < 0.1 &&
			Perlin::peaks.GetPercentBeneathThreshold(x, y) > 2 * Perlin::oceans.GetPercentAboveThreshold(x, y))
		{
			if (!Perlin::randomLakes.UnderThreshold(x, y))
			{
				SetOcean();
				return;
			}
			if (!Perlin::randomHills.UnderThreshold(x, y))
			{
				SetOcean();
				return;
			}
			if (Perlin::randomPasses.UnderThreshold(x, y) && Perlin::randomPasses.GetPercentBeneathThreshold(x, y) > 0.1)
			{
				SetHills();
				return;
			}
			SetPeaks();
			return;
		}
		if (Perlin::mountains.UnderThreshold(x, y) &&
			Perlin::mountains.GetPercentBeneathThreshold(x, y) > 0.1 &&
			Perlin::oceans.GetPercentAboveThreshold(x, y) < 0.15 &&
			Perlin::mountains.GetPercentBeneathThreshold(x, y) > 2 * Perlin::oceans.GetPercentAboveThreshold(x, y))
		{
			if (!Perlin::randomLakes.UnderThreshold(x, y))
			{
				SetOcean();
				return;
			}
			if (!Perlin::randomHills.UnderThreshold(x, y))
			{
				SetOcean();
				return;
			}
			if (Perlin::randomPasses.UnderThreshold(x, y))
			{
				SetHills();
				return;
			}
			SetMountains();
			return;
		}
		SetOcean();
		return;
	}
	if (Perlin::peaks.UnderThreshold(x, y))
	{
		if (!Perlin::randomPasses.UnderThreshold(x, y))
		{
			SetPeaks();
			return;
		}
		else
		{
			SetHills();
			return;
		}
	}
	if (Perlin::mountains.UnderThreshold(x, y))
	{
		if (!Perlin::randomPasses.UnderThreshold(x, y))
		{
			SetMountains();
			return;
		}
		else
		{
			SetHills();
			return;
		}
	}
	if (Perlin::foothills.UnderThreshold(x, y) && !Perlin::randomPasses.UnderThreshold(x, y))
	{
		SetHills();
		return;
	}
	if (!Perlin::randomLakes.UnderThreshold(x, y))
	{
		SetLakes();
		return;
	}

	if (!Perlin::randomHills.UnderThreshold(x, y)
		&& Perlin::oceans.GetPercentBeneathThreshold(x, y) > 0.05
		&& Perlin::randomLakes.GetPercentBeneathThreshold(x, y) > 0.08)
	{
		SetHills();
		return;
	}

	SetFlats();
	return;
}

double 
SamplePoint::GetBaseSedimentDepth()
{
	if (type.IsTerrainOfType(TerrainTemplate::PEAKS))
		return 0;
	else if (type.IsTerrainOfType(TerrainTemplate::MOUNTAINS))
		return 3;
	else if (type.IsTerrainOfType(TerrainTemplate::ROUGH))
		return 12;
	else
		return 60;
}

vector<double> 
SamplePoint::GetPerlinElevDiffs()
{
	if (IsWaterPoint())
		return vector<double>{0, 0};
	else if (type.IsTerrainOfType(TerrainTemplate::PEAKS))
		return vector<double>{0, 1};
	else if (type.IsTerrainOfType(TerrainTemplate::MOUNTAINS))
		return vector<double>{0.1, 0.9};
	else if (type.IsTerrainOfType(TerrainTemplate::HILLS))
		return vector<double>{0.4, 0.6};
	else
		return vector<double>{1, 0};
}

glm::vec3 
SamplePoint::GetTerrainBasedColor()
{
	if(type.IsTerrainOfType(TerrainTemplate::OCEAN))
		return glm::vec3(14./255, 21./255, 110./255);
	if (type.IsTerrainOfType(TerrainTemplate::LAKE) && !type.IsTerrainOfType(TerrainTemplate::FLAT))
		return glm::vec3(71./255, 94./255, 99./255);
	if (type.IsTerrainOfType(TerrainTemplate::PEAKS))
		return glm::vec3(171./255, 156./255, 135./255);
	if (type.IsTerrainOfType(TerrainTemplate::MOUNTAINS))
		return glm::vec3(148./255, 126./255, 40./255);
	if (type.IsTerrainOfType(TerrainTemplate::HILLS))
		return glm::vec3(156./255, 158./255, 93./255);
	if (type.IsTerrainOfType(TerrainTemplate::LAKE))
		return glm::vec3(18./255, 146./255, 201./255);
	if (type.IsTerrainOfType(TerrainTemplate::FLAT))
		return glm::vec3(124./255, 166./255, 88./255);
	return glm::vec3(1.0, 0.0, 0.0);
}

bool
SamplePoint::IsInlandLake() {
	if (IsOcean())
		return false;
	return type.IsTerrainOfType(TerrainTemplate::LAKE);
}
double 
SamplePoint::GetTectonicUplift()
{
	if (tectonicUpliftOverride != -1)
		return tectonicUpliftOverride;
	//Paper uses 5.01 * 10^-4 (ie 0.0001)
	double adjust = Perlin::upliftAdjust.Get(x, y);
	if (type.IsTerrainOfType(TerrainTemplate::PEAKS))
		return (6.01 + 1.5 * adjust) * 0.0001;
	else if (type.IsTerrainOfType(TerrainTemplate::MOUNTAINS))
		return (4.52 + 1.5 * adjust) * 0.0001;
	else if (type.IsTerrainOfType(TerrainTemplate::HILLS))
		return (adjust > 0 ? (1.25 + 2.5 * adjust) : (1.25 + 0.75 * adjust)) * 0.0001;
	else if (IsOcean())
		return 0;
	else
		return (4.15 + 3.5 * adjust) * 0.00001;
}