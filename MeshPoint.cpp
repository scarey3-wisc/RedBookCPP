#include "MeshPoint.h"
#include <algorithm>

/*public MeshPoint(DataInputStream dis)
{
	try
	{
		this.x = dis.readDouble();
		this.y = dis.readDouble();
		double elev = dis.readDouble();
		SetElevation(elev);
	}
	catch (IOException e)
	{
		e.printStackTrace();
	}
	adjacent = new ConcurrentHashMap<MeshPoint, MeshConnection>();
	searchID = 0;
	ResetContainerIndex();
	ResetDrainage();
	ResetTag();
}
public MeshPoint(Iterator<String> desc)
{
	this.x = Double.parseDouble(desc.next());
	this.y = Double.parseDouble(desc.next());
	double elev = Double.parseDouble(desc.next());
	SetElevation(elev);
	adjacent = new ConcurrentHashMap<MeshPoint, MeshConnection>();
	searchID = 0;
	ResetContainerIndex();
	ResetDrainage();
	ResetTag();
}
public static boolean ConsumeDescription(DataInputStream dis)
{
	try {
		dis.readDouble();
		dis.readDouble();
		dis.readDouble();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public boolean WriteDescription(DataOutputStream dos)
{
	try {
		dos.writeDouble(x);
		dos.writeDouble(y);
		dos.writeDouble(currElev);
		return true;
	}
	catch (IOException e) {
		return false;
	}
}

public String GetDescription()
{
	String desc = "";
	desc += Double.toString(x) + " ";
	desc += Double.toString(y) + " ";
	desc += Double.toString(currElev);
	return desc;
}*/

void 
MeshPoint::InitDirectlyAdjacent()
{
	ResetDirectlyAdjacent();
	for (auto p : adjacent)
	{

		MeshPoint* adj = GetClosestDirectAdjacency(p.first);
		directlyAdjacent.push_back(adj);
	}
}

void 
MeshPoint::CalculateDrainageArea()
{
	totalDrainageArea = personalDrainageArea;
	for (MeshPoint* mp : sources)
	{
		totalDrainageArea += mp->GetDrainageArea();
	}
}



void 
MeshPoint::AssignDrainage()
{
	for (MeshPoint* adj : directlyAdjacent)
	{
		if (adj->GetElevation() >= GetElevation())
			continue;
		if (sink == nullptr)
		{
			sink = adj;
			continue;
		}
		if (IsWaterPoint())
		{
			if (sink->IsWaterPoint() && !adj->IsWaterPoint())
				continue;
			if (adj->IsWaterPoint() && !sink->IsWaterPoint())
			{
				sink = adj;
				continue;
			}
		}
		if (adj->GetElevation() < sink->GetElevation())
			sink = adj;
	}
	if (sink != nullptr)
		sink->sources.push_back(this);
}

//the forced drainage doesn't need to be adjacent to us
//this is primarily in the case of "local minima", where
//we need to teleport water from the minima to the nearest
//pass to keep water running downhill. 
void 
MeshPoint::ForceAssignDrainage(MeshPoint* targetSink)
{
	sink = targetSink;
	if (sink != nullptr)
		sink->sources.push_back(this);
}

void 
MeshPoint::ResetDrainage()
{
	sink = nullptr;
	sources.clear();
	totalDrainageArea = 0;
}

//Just looking at my direct adjacencies, what drains into me?
void 
MeshPoint::CalculatePersonalDrainageArea()
{
	personalDrainageArea = 0;
	std::vector<MeshPoint*> directlyAdjacent;
	for (auto p : adjacent)
	{
		MeshPoint* dir = GetClosestDirectAdjacency(p.first);
		if (dir != nullptr)
			directlyAdjacent.push_back(dir);
	}
	std::sort(directlyAdjacent.begin(), directlyAdjacent.end(),
		[this](const MeshPoint*& o1, const MeshPoint*& o2) {
			double theta1 = std::atan2(o1->x - x, o1->y - y);
			double theta2 = std::atan2(o2->x - x, o2->y - y);
			if (theta1 < theta2)
				return -1;
			if (theta1 > theta2)
				return 1;
			return 0;
		});
	for (int i = 0; i < directlyAdjacent.size(); i++)
	{
		MeshPoint* one = directlyAdjacent[i];
		MeshPoint* two = directlyAdjacent[(i + 1) % directlyAdjacent.size()];
		double area = VoronoiAlgorithms.TriangleArea(this, one, two);
		personalDrainageArea += area / 3;
	}
}

bool 
MeshPoint::HasZeroElevDiffs()
{
	std::vector<double> perlinElevDiffs = GetPerlinElevDiffs();
	for (int i = 0; i < perlinElevDiffs.size(); i++)
	{
		if (perlinElevDiffs[i] != 0)
			return false;
	}
	return true;
}

void 
MeshPoint::ResetAdjacencies()
{
	for (auto e : adjacent)
	{
		if (e.second->MidInitialized())
			e.second->GetMid().ResetAdjacencies();
		e.first->adjacent.erase(this);
	}
	adjacent.clear();
}

void 
MeshPoint::ForceRemoveAdjacency(MeshPoint* target)
{
	adjacent.erase(target);
}

void 
MeshPoint::ForceOneWayAdjacency(MeshPoint* target)
{
	MeshConnection* con = new MeshConnection(this, target);
	adjacent[target] = con;
}

void 
MeshPoint::MarkAdjacent(MeshPoint* p)
{
	if (adjacent.find(p) != adjacent.end() && p->adjacent.find(this) != p->adjacent.end())
		return;
	if (p->adjacent.find(this) != p->adjacent.end())
		p->adjacent.erase(this);
	if (adjacent.find(p) != adjacent.end())
		adjacent.erase(p);
	MeshConnection* con = new MeshConnection(this, p);
	adjacent[p] = con;
	p->adjacent[this] = con;
}

MeshConnection*
MeshPoint::GetConnection(MeshPoint* adj)
{
	if (adjacent.find(adj) == adjacent.end())
		return nullptr;
	return adjacent[adj];
}

//The idea here is that our adjacency list might
//include points that we aren't "directly adjacent to"
//because the connection has a midpoint; we're directly
//adjacent to that midpoint. This descends the tree to
//those midpoints. 
MeshPoint* 
MeshPoint::GetClosestDirectAdjacency(MeshPoint* adj)
{
	MeshConnection con = GetConnection(adj);
	if (con.GetMid() == null)
		return adj;
	while (con.GetMid() != null)
	{
		adj = con.GetMid();
		con = adj.GetConnection(this);
	}
	return adj;
}

double 
MeshPoint::ConvertVoronoiDistToMeters(double d)
{
	return d * RegionalMap.DIMENSION * LocalMap.METER_DIM;
}