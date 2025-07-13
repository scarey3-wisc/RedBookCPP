#include "VoronoiAlgorithms.h"
#include "MeshConnection.h"
#include "Perlin.h"
#include <queue>
#include <limits>
#include "MathToolkit.h"


using namespace std;

void 
VoronoiAlgorithms::EnsureRiversFlowDownhill(vector<SamplePoint*>& continent)
{
	for (SamplePoint* s : continent)
	{
		if (s->GetRiverOutlet() == nullptr)
			continue;
		MeshConnection* mc = s->GetConnection(s->GetRiverOutlet());
		if (mc == nullptr)
			continue;
		if (!mc->MidInitialized())
			continue;
		if (s->GetRiverOutlet()->IsOcean() || s->GetRiverOutlet()->type.IsTerrainOfType(TerrainTemplate::LAKE))
			continue;

		MeshMidPoint* mid = mc->GetMid();
		mid->SetElevation((s->GetElevation() + s->GetRiverOutlet()->GetElevation()) / 2);
		double elevChange = 0;
		for (int i = 0; i < Perlin::elevDeltaLength; i++)
		{
			double del = abs(Perlin::elevDeltas[i].Get(mid->x, mid->y));
			double ctr = mid->GetPerlinElevDiffs()[i];
			double amp = Perlin::elevDeltaScales[i];
			double elevIncrease = amp * ctr * del;
			elevChange += elevIncrease;
		}
		mid->SetElevation(mid->GetElevation() - elevChange);
	}
	for (SamplePoint* s : continent)
	{
		if (s->GetRiverInlets().size() != 0)
		{
			double elevChange = 0;
			for (int i = 0; i < Perlin::elevDeltaLength; i++)
			{
				double del = abs(Perlin::elevDeltas[i].Get(s->x, s->y));
				double ctr = s->GetPerlinElevDiffs()[i];
				double amp = Perlin::elevDeltaScales[i];
				double elevIncrease = amp * ctr * del;
				elevChange += elevIncrease;
			}
			s->SetElevation(s->GetElevation() - elevChange);
		}
	}
}
void 
VoronoiAlgorithms::IncreaseFractureLevel(vector<SamplePoint*>& area)
{
	vector<vector<MeshPoint*>> allDetailLevels;
	vector<MeshPoint*> currentDetailLevel;
	for (SamplePoint* s : area)
		currentDetailLevel.push_back(s);

	vector<MeshPoint*> nextDetailLevel;
	while (true)
	{
		for (MeshPoint* mp : currentDetailLevel)
		{
			mp->ForEachAdjacent([&nextDetailLevel, mp](MeshPoint* adj) {
				MeshConnection* con = mp->GetConnection(adj);
				if (con == nullptr)
					return;
				if (con->GetMid() != nullptr)
					nextDetailLevel.push_back(con->GetMid());
			});
		}
		allDetailLevels.push_back(currentDetailLevel);
		currentDetailLevel = nextDetailLevel;
		nextDetailLevel.clear();
		if (currentDetailLevel.size() == 0)
			break;
	}
	for (auto detLevel : allDetailLevels)
	{
		for (MeshPoint* mp : detLevel)
		{
			mp->ForEachAdjacent([mp](MeshPoint* adj) {
				MeshConnection* con = mp->GetConnection(adj);
				if (con == nullptr)
					return;
				if (con->GetMid() != nullptr)
					return;

				MeshMidPoint* newMid = con->GetMid(true);
				if (mp->IsWaterPoint() && adj->IsWaterPoint())
				{
					if (mp->IsOcean())
					{
						newMid->SetOcean();
						newMid->SetElevation(mp->GetElevation());
						newMid->CopyPerlinElevDiffs(mp->GetPerlinElevDiffs());
					}
					else if (adj->IsOcean())
					{
						newMid->SetOcean();
						newMid->SetElevation(adj->GetElevation());
						newMid->CopyPerlinElevDiffs(adj->GetPerlinElevDiffs());
					}
					else
						newMid->SetInlandLake();
				}
				else if (mp->IsWaterPoint())
				{
					newMid->SetElevation(mp->GetElevation());
					newMid->CopyPerlinElevDiffs(mp->GetPerlinElevDiffs());
					if (mp->IsOcean())
						newMid->SetOcean();
					else
						newMid->SetInlandLake();
				}
				else if (adj->IsWaterPoint())
				{
					newMid->SetElevation(adj->GetElevation());
					newMid->CopyPerlinElevDiffs(adj->GetPerlinElevDiffs());
					if (adj->IsOcean())
						newMid->SetOcean();
					else
						newMid->SetInlandLake();
				}
			});
		}
	}
}
void 
VoronoiAlgorithms::SetRiverFlow(vector<SamplePoint*>& continent)
{
	VoronoiAlgorithms.SortByElevation(continent);
	for (SamplePoint* c : continent)
	{
		c->ResetRiver();
	}
	for (int i = continent.size() - 1; i >= 0; i--)
	{
		SamplePoint* p = continent[i];
		if (p->RiverFlowProcessed())
			continue;
		if (i == 0)
		{
			continent[i]->SendRiverFlow();
			break;
		}

		//We're checking for a lake block; a block of points all at the same elev
		//j is the index of the next point that goes lower
		bool atOcean = false;
		int j = i - 1;
		while (continent[j]->GetElevation() == p->GetElevation())
		{
			j--;
			if (j < 0)
			{
				atOcean = true;
				j++;
				break;
			}
		}
		if (atOcean)
			continue;

		if (j == i - 1)
		{
			continent[i]->SendRiverFlow();
			continue;
		}

		//These are locations on our lake block which have a way down
		vector<SamplePoint*> found;
		queue<SamplePoint*> horizon;
		for (int k = i; k != j; k--)
		{
			if (continent[k]->GetWayDownhill(true, false) != nullptr)
				horizon.push(continent[k]);
		}
		while (!horizon.empty())
		{
			SamplePoint* curr = horizon.front();
			horizon.pop();
			found.push_back(curr);
			for (SamplePoint* adj : curr->GrabAdjacentUnassignedRiverFlows(false))
				horizon.push(adj);
		}
		for (int k = found.size() - 1; k >= 0; k--)
		{
			found[k]->SendFlowToOutlet();
		}
	}
}
void 
VoronoiAlgorithms::ConvertCoastalLakeToOcean(vector<SamplePoint*>& continentArea)
{
	vector<SamplePoint*> lakeShore = VoronoiAlgorithms::FindBoundaryPoints(continentArea, TerrainTemplate::LAKE);
	vector<vector<SamplePoint*>> lakes =
		VoronoiAlgorithms::FindTypeClumps(lakeShore, TerrainTemplate::LAKE, numeric_limits<int>::max());
	for (vector<SamplePoint*> lake : lakes)
	{
		bool isCoastalLake = false;
		for (SamplePoint* l : lake)
		{
			for (SamplePoint* a : l->GetAdjacentSamples())
			{
				if (a->type.IsTerrainOfType(TerrainTemplate::OCEAN))
					isCoastalLake = true;
			}
		}
		if (isCoastalLake)
			for (SamplePoint* l : lake)
				l->MakeOcean();
	}
}
void 
VoronoiAlgorithms::ConvertSeasToLakes(vector<SamplePoint*>& continentCoast, int maxLakeSize)
{
	vector<vector<SamplePoint*>> seas =
		VoronoiAlgorithms.FindTypeClumps(continentCoast, TerrainTemplate::OCEAN, maxLakeSize);
	for (vector<SamplePoint*> gp : seas)
	{
		for (SamplePoint* v : gp)
		{
			v->MakeLake();
		}
	}
}

glm::dvec3 
VoronoiAlgorithms::BarycentricCoordinatesDY(double x, double y, array<MeshPoint*, 3> vertices)
{
	MeshPoint* a = vertices[0];
	MeshPoint* b = vertices[1];
	MeshPoint* c = vertices[2];
	double abc = MathToolkit::SignedTriangleArea(a->x, a->y, b->x, b->y, c->x, c->y);
	double pbc = SignedTriangleAreaDY(x, y, b->x, b->y, c->x, c->y, 1);
	double apc = SignedTriangleAreaDY(a->x, a->y, x, y, c->x, c->y, 2);
	double abp = SignedTriangleAreaDY(a->x, a->y, b->x, b->y, x, y, 3);
	return glm::dvec3(pbc / abc, apc / abc, abp / abc);
}
double 
VoronoiAlgorithms::SignedTriangleAreaDY(double x1, double y1, double x2, double y2, double x3, double y3, int whichCoord)
{
	if (whichCoord == 1)
	{
		return 0.5 * (x3 - x2);
	}
	if (whichCoord == 2)
	{
		return 0.5 * (x1 - x3);
	}
	if (whichCoord == 3)
	{
		return 0.5 * (x2 - x1);
	}
	return 0;
}
glm::dvec3
VoronoiAlgorithms::BarycentricCoordinatesDX(double x, double y, array<MeshPoint*, 3> vertices)
{
	MeshPoint* a = vertices[0];
	MeshPoint* b = vertices[1];
	MeshPoint* c = vertices[2];
	double abc = MathToolkit::SignedTriangleArea(a->x, a->y, b->x, b->y, c->x, c->y);
	double pbc = SignedTriangleAreaDX(x, y, b->x, b->y, c->x, c->y, 1);
	double apc = SignedTriangleAreaDX(a->x, a->y, x, y, c->x, c->y, 2);
	double abp = SignedTriangleAreaDX(a->x, a->y, b->x, b->y, x, y, 3);
	return glm::dvec3(pbc / abc, apc / abc, abp / abc);
}
double 
VoronoiAlgorithms::SignedTriangleAreaDX(double x1, double y1, double x2, double y2, double x3, double y3, int whichCoord)
{
	if (whichCoord == 1)
	{
		return 0.5 * (y2 - y3);
	}
	if (whichCoord == 2)
	{
		return 0.5 * (y3 - y1);
	}
	if (whichCoord == 3)
	{
		return 0.5 * (y1 - y2);
	}
	return 0;
}

glm::dvec3
VoronoiAlgorithms::BarycentricCoordinates(double x, double y, array<SamplePoint*, 3> vertices)
{
	return BarycentricCoordinates(x, y, array<MeshPoint*, 3>{vertices[0], vertices[1], vertices[2]});
}

glm::dvec3 
VoronoiAlgorithms::BarycentricCoordinates(double x, double y, array<MeshPoint*, 3> vertices)
{
	MeshPoint* a = vertices[0];
	MeshPoint* b = vertices[1];
	MeshPoint* c = vertices[2];
	double abc = MathToolkit::SignedTriangleArea(a->x, a->y, b->x, b->y, c->x, c->y);
	double pbc = MathToolkit::SignedTriangleArea(x, y, b->x, b->y, c->x, c->y);
	double apc = MathToolkit::SignedTriangleArea(a->x, a->y, x, y, c->x, c->y);
	double abp = MathToolkit::SignedTriangleArea(a->x, a->y, b->x, b->y, x, y);
	return glm::dvec3(pbc / abc, apc / abc, abp / abc);
}
double 
VoronoiAlgorithms::TriangleArea(MeshPoint* a, MeshPoint* b, MeshPoint* c)
{
	double signedArea = MathToolkit::SignedTriangleArea(
		a->x, a->y,
		b->x, b->y,
		c->x, c->y);
	return abs(signedArea);
}

double 
VoronoiAlgorithms::TriangleArea(array<MeshPoint*, 3> vertices)
{
	double signedArea = MathToolkit::SignedTriangleArea(
		vertices[0]->x, vertices[0]->y,
		vertices[1]->x, vertices[1]->y,
		vertices[2]->x, vertices[2]->y);
	return abs(signedArea);
}

array<MeshPoint*, 3>
VoronoiAlgorithms::FindContainingTriangle(double x, double y, SamplePoint* seed)
{
	array<SamplePoint*, 3> start = FindContainingSampleTriangle(x, y, seed);
	array<MeshPoint*, 3> curr { start[0], start[1], start[2] };
	if (curr[0] == nullptr || curr[1] == nullptr || curr[2] == nullptr)
		return array<MeshPoint*, 3> {nullptr, nullptr, nullptr};
	while (curr[0] != nullptr)
	{
		bool subtriangleFound = false;
		glm::dvec2 coords = BarycentricCoordinates(x, y, curr);
		for (int i = 0; i < 3; i++)
		{
			if (coords[i] >= 0.5)
			{
				MeshConnection* a = MeshConnection::FindConnection(curr[i], curr[(i + 1) % 3]);
				MeshConnection* b = MeshConnection::FindConnection(curr[i], curr[(i + 2) % 3]);
				if (a == nullptr || b == nullptr)
					return curr;
				if (!a->MidInitialized() && !b->MidInitialized())
					return curr;
				else if (a->MidInitialized() && b->MidInitialized())
				{
					curr = array<MeshPoint*, 3>{ curr[i], a->GetMid(), b->GetMid() };
					subtriangleFound = true;
				}
				else if (a->MidInitialized() && !b->MidInitialized())
				{
					MeshPoint* bMid = b->GetMidOrTempMid();
					curr = array<MeshPoint*, 3>{ curr[i], a->GetMid(), bMid };
					subtriangleFound = true;
				}
				else if (!a->MidInitialized() && b->MidInitialized())
				{
					MeshPoint* aMid = a->GetMidOrTempMid();
					curr = array<MeshPoint*, 3>{ curr[i], aMid, b->GetMid() };
					subtriangleFound = true;
				}
			}
		}
		if (!subtriangleFound)
		{
			MeshConnection* a = MeshConnection::FindConnection(curr[0], curr[1]);
			MeshConnection* b = MeshConnection::FindConnection(curr[0], curr[2]);
			MeshConnection* c = MeshConnection::FindConnection(curr[1], curr[2]);
			if (a == nullptr || b == nullptr || c == nullptr)
				return curr;
			int numMids = 0;
			if (a->MidInitialized())
				numMids++;
			if (b->MidInitialized())
				numMids++;
			if (c->MidInitialized())
				numMids++;
			if (numMids == 3)
				curr = array<MeshPoint*, 3>{ a->GetMid(), b->GetMid(), c->GetMid() };
			else if (numMids == 2)
			{
				MeshPoint* aMid = a->GetMid();
				if (aMid == nullptr)
					aMid = new MeshMidPoint(a, false);
				MeshPoint* bMid = b->GetMid();
				if (bMid == nullptr)
					bMid = new MeshMidPoint(b, false);
				MeshPoint* cMid = c->GetMid();
				if (cMid == nullptr)
					cMid = new MeshMidPoint(c, false);
				curr = array<MeshPoint*, 3>{ aMid, bMid, cMid };
			}
			else
				return curr;
		}
	}
	return curr;

}
array<SamplePoint*, 3>
VoronoiAlgorithms::FindContainingSampleTriangle(double x, double y, SamplePoint* seed)
{
	array<SamplePoint*, 3> results;
	results[0] = seed;
	int attempts = 0;
	while (attempts < 5)
	{
		attempts++;
		SamplePoint* curr = results[0];
		double vx = x - curr->x;
		double vy = y - curr->y;
		double vD = sqrt(vx * vx + vy * vy);
		vx /= vD;
		vy /= vD;
		SamplePoint* best = nullptr;
		double bestDot = -1;
		for (SamplePoint* a : curr->GetAdjacentSamples())
		{
			double avx = a->x - curr->x;
			double avy = a->y - curr->y;
			double avD = sqrt(avx * avx + avy * avy);
			avx /= avD;
			avy /= avD;
			double dot = vx * avx + vy * avy;
			if (best == nullptr || dot > bestDot)
			{
				best = a;
				bestDot = dot;
			}
		}
		SamplePoint* bestOtherSide = nullptr;
		double bestOtherDot = -1;
		for (SamplePoint* b : curr->GetAdjacentSamples())
		{
			if (b == nullptr)
				continue;
			if (b == best)
				continue;
			double bvx = b->x - curr->x;
			double bvy = b->y - curr->y;
			double bvD = sqrt(bvx * bvx + bvy * bvy);
			bvx /= bvD;
			bvy /= bvD;
			double dot = vx * bvx + vy * bvy;

			double avx = best->x - curr->x;
			double avy = best->y - curr->y;
			double avD = sqrt(avx * avx + avy * avy);
			avx /= avD;
			avy /= avD;

			//using the cross product to check that the vector to our point lies between the two edge vectors
			double aXv = avy * vx - avx * vy;
			double aXb = avy * bvx - avx * bvy;
			double bXv = bvy * vx - bvx * vy;
			if (aXv * aXb >= 0 && bXv * aXb * -1 >= 0)
			{
				if (bestOtherSide == nullptr || dot > bestOtherDot)
				{
					bestOtherSide = b;
					bestOtherDot = dot;
				}
			}
		}
		if (best == nullptr || bestOtherSide == nullptr)
		{
			results[0] = nullptr;
			return results;
		}
		if (curr->DistTo(best) < curr->DistTo(bestOtherSide))
		{
			results[1] = best;
			results[2] = bestOtherSide;
		}
		else
		{
			results[1] = bestOtherSide;
			results[2] = best;
		}
		//so now we have a candidate triangle: are we inside it?
		glm::dvec3 bcCoords = BarycentricCoordinates(x, y, results);
		if (bcCoords[0] < 0 || bcCoords[1] < 0 || bcCoords[2] < 0)
		{
			results[0] = results[1];
			results[1] = nullptr;
			results[2] = nullptr;
			//return null;
		}
		else
		{
			return results;
		}
	}
	results[0] = nullptr;
	results[1] = nullptr;
	results[2] = nullptr;
	return results;
}
private static void ClearHeights(ArrayList<SamplePoint> coastalPoints)
{
	LinkedList<SamplePoint> frontier = new LinkedList<SamplePoint>();
	SamplePoint.StartNewSearch();
	for (SamplePoint vp : coastalPoints)
	{
		if (vp.IsOcean() || vp.Reached())
			continue;
		vp.MarkAsReached();
		frontier.addLast(vp);
	}
	while (!frontier.isEmpty())
	{
		SamplePoint curr = frontier.removeFirst();
		curr.ResetElevation();
		for (SamplePoint vp : curr.GetAdjacentSamples())
		{
			if (vp.IsOcean() || vp.Reached())
				continue;
			frontier.addLast(vp);
			vp.MarkAsReached();
		}
	}
}
public static void AssignHeights(ArrayList<SamplePoint> coastalPoints)
{
	ClearHeights(coastalPoints);
	PriorityQueue<VoronoiNode> frontier = new PriorityQueue<VoronoiNode>();
	SamplePoint.StartNewSearch();
	for (SamplePoint vp : coastalPoints)
	{
		if (vp.IsOcean())
			continue;
		double distToOcean = -1;
		for (SamplePoint adj : vp.GetAdjacentSamples())
		{
			if (adj.IsOcean())
			{
				double dist = vp.DistTo(adj);
				if (distToOcean == -1 || dist < distToOcean)
					distToOcean = dist;
			}
		}
		if (distToOcean == -1)
			continue;
		distToOcean = SamplePoint.ConvertVoronoiDistToMeters(distToOcean / 2);
		double grade = QueryGrade(TerrainType.OCEAN, vp);

		double delta = distToOcean * grade;
		vp.SetElevation(delta);
		frontier.offer(new VoronoiNode(vp, vp.GetElevation()));
	}
	while (!frontier.isEmpty())
	{
		VoronoiNode curr = frontier.poll();
		if (curr.v.Reached())
			continue;
		curr.v.MarkAsReached();
		for (SamplePoint adj : curr.v.GetAdjacentSamples())
		{
			if (adj.IsOcean() || adj.Reached())
				continue;
			double distTo = SamplePoint.ConvertVoronoiDistToMeters(curr.v.DistTo(adj));
			if (curr.v.type.IsTerrainOfType(TerrainType.LAKE) || curr.v.type.IsTerrainOfType(TerrainTemplate.OCEAN))
				distTo /= 2;

			double grade = QueryGrade(curr.v.type, adj);

			double candDelta = distTo * grade;
			double candValue = curr.v.GetElevation() + candDelta;
			if (candValue < adj.GetElevation())
			{
				adj.SetElevation(candValue);
				frontier.offer(new VoronoiNode(adj, adj.GetElevation()));
			}

		}
	}
}
public static void SortByElevation(ArrayList< ? extends MeshPoint> continent)
{
	Comparator<MeshPoint> compare = new Comparator<MeshPoint>()
	{

		@Override
			public int compare(MeshPoint o1, MeshPoint o2) {
			if (o1.GetElevation() < o2.GetElevation())
				return -1;
			if (o1.GetElevation() > o2.GetElevation())
				return 1;
			return 0;
		}
	};
	continent.sort(compare);
}

public static ArrayList<SamplePoint> FindAllOfType(ArrayList<SamplePoint> area, TerrainTemplate target)
{
	ArrayList<SamplePoint> found = new ArrayList<SamplePoint>();
	for (SamplePoint sp : area)
		if (sp.type.IsTerrainOfType(target))
			found.add(sp);
	return found;
}
//NOTE: the seeds are assumed to be *adjacent* to Terrain of the target type.
//The classic use case here is having coastal 
public static ArrayList<ArrayList<SamplePoint>> FindTypeClumps(ArrayList<SamplePoint> seeds, TerrainTemplate target, int maxNum)
{
	ArrayList<ArrayList<SamplePoint>> found = new ArrayList<ArrayList<SamplePoint>>();
	SamplePoint.StartNewSearch();
	for (SamplePoint cs : seeds)
	{
		for (SamplePoint ts : cs.GetAdjacentSamples())
		{
			ArrayList<SamplePoint> group = new ArrayList<SamplePoint>();
			if (ts.Reached())
				continue;
			if (!ts.type.IsTerrainOfType(target))
				continue;
			LinkedList<SamplePoint> horizon = new LinkedList<SamplePoint>();
			ts.MarkAsReached();
			horizon.addLast(ts);
			while (!horizon.isEmpty())
			{
				SamplePoint curr = horizon.removeFirst();
				group.add(curr);
				for (SamplePoint adj : curr.GetAdjacentSamples())
				{
					if (adj.Reached())
						continue;
					if (!adj.type.IsTerrainOfType(target))
						continue;
					adj.MarkAsReached();
					horizon.addLast(adj);
				}
			}
			if (group.size() <= maxNum)
				found.add(group);
		}
	}
	return found;
}
public static ArrayList<SamplePoint> FindAllWithinBoundary(SamplePoint seed, TerrainTemplate boundaryType)
{
	ArrayList<SamplePoint> listedSeeds = new ArrayList<SamplePoint>();
	listedSeeds.add(seed);
	return FindAllWithinBoundary(listedSeeds, boundaryType);
}
public static ArrayList<SamplePoint> FindAllWithinBoundary(ArrayList<SamplePoint> seeds, TerrainTemplate boundaryType)
{
	SamplePoint.StartNewSearch();
	ArrayList<SamplePoint> found = new ArrayList<SamplePoint>();
	LinkedList<SamplePoint> horizon = new LinkedList<SamplePoint>();
	for (SamplePoint seed : seeds)
	{
		if (seed.type.IsTerrainOfType(boundaryType))
			continue;
		if (seed.Reached())
			continue;
		horizon.add(seed);
		seed.MarkAsReached();
	}

	while (!horizon.isEmpty())
	{
		SamplePoint curr = horizon.removeFirst();
		for (SamplePoint vp : curr.GetAdjacentSamples())
		{
			if (vp.Reached())
				continue;
			if (!vp.type.IsTerrainOfType(boundaryType))
			{
				vp.MarkAsReached();
				horizon.addLast(vp);
			}
		}
		found.add(curr);
	}
	return found;
}
public static ArrayList<SamplePoint> FindBoundaryPoints(SamplePoint seed, TerrainTemplate boundaryTo)
{
	ArrayList<SamplePoint> listedSeeds = new ArrayList<SamplePoint>();
	listedSeeds.add(seed);
	return FindBoundaryPoints(listedSeeds, boundaryTo);
}
public static ArrayList<SamplePoint> FindBoundaryPoints(ArrayList<SamplePoint> seeds, TerrainTemplate boundaryTo)
{
	SamplePoint.StartNewSearch();
	ArrayList<SamplePoint> found = new ArrayList<SamplePoint>();
	LinkedList<SamplePoint> horizon = new LinkedList<SamplePoint>();
	for (SamplePoint seed : seeds)
	{
		if (!seed.Reached())
			horizon.add(seed);
		seed.MarkAsReached();
	}
	while (!horizon.isEmpty())
	{
		SamplePoint curr = horizon.removeFirst();
		if (curr.type.IsTerrainOfType(boundaryTo))
		{
			for (SamplePoint vp : curr.GetAdjacentSamples())
			{
				if (vp.Reached())
					continue;
				vp.MarkAsReached();
				if (vp.type.IsTerrainOfType(boundaryTo))
					horizon.addLast(vp);
				else
					found.add(vp);
			}
		}
		else
		{
			boolean amBoundary = false;
			for (SamplePoint vp : curr.GetAdjacentSamples())
			{
				if (vp.type.IsTerrainOfType(boundaryTo))
					amBoundary = true;
				if (vp.Reached())
					continue;
				if (!vp.type.IsTerrainOfType(boundaryTo))
				{
					vp.MarkAsReached();
					horizon.addLast(vp);
				}
			}
			if (amBoundary)
				found.add(curr);
		}
	}
	return found;
}
private static class VoronoiNode implements Comparable<VoronoiNode>
{
	private double dist;
	public SamplePoint v;
	public VoronoiNode(SamplePoint vp, double dist)
	{
		v = vp;
		this.dist = dist;
	}
	@Override
		public int compareTo(VoronoiNode o) {
		if (dist > o.dist)
			return 1;
		if (dist < o.dist)
			return -1;
		return 0;
	}
}
private static double QueryGrade(TerrainTemplate from, SamplePoint to)
{
	double minGrade = QueryGrade(from, to.type, true);
	double maxGrade = QueryGrade(from, to.type, false);
	double perlinPush = Perlin.minMaxSelector.Get(to.x, to.y);
	perlinPush += 1;
	perlinPush /= 2;
	if (perlinPush < 0)
		perlinPush = 0;
	if (perlinPush > 1)
		perlinPush = 1;
	double grade = perlinPush * maxGrade + (1 - perlinPush) * minGrade;
	return grade;
}
private static double QueryGrade(TerrainTemplate from, TerrainTemplate to, boolean min)
{
	int fromIndex = GradeQueryIndex(from);
	int toIndex = GradeQueryIndex(to);
	if (min)
		return GRADE_MIN_MATRIX[fromIndex][toIndex];
	else
		return GRADE_MAX_MATRIX[fromIndex][toIndex];
}
private static int GradeQueryIndex(TerrainTemplate type)
{
	if (type.IsTerrainOfType(TerrainTemplate.PEAKS))
		return 4;
	else if (type.IsTerrainOfType(TerrainTemplate.MOUNTAINS))
		return 3;
	else if (type.IsTerrainOfType(TerrainTemplate.HILLS))
		return 2;
	else if (type.IsTerrainOfType(TerrainTemplate.OCEAN))
		return 1;
	else if (type.IsTerrainOfType(TerrainTemplate.LAKE))
		return 1;
	else
		return 0;
}
private static double[][] GRADE_MAX_MATRIX = new double[][]
{
	{0.0040, 0.0000, 0.0250, 0.1800, 0.3000}, //From flatland to...
	{0.0040, 0.0000, 0.1500, 0.5000, 0.5000}, //From ocean to...
	{0.0005, 0.0000, 0.0200, 0.1800, 0.3000}, //From hills to...
	{0.0005, 0.0000, 0.0100, 0.0300, 0.3000}, //From mountains to...
	{0.0005, 0.0000, 0.0100, 0.0100, 0.3000}  //From peaks to...
};
private static double[][] GRADE_MIN_MATRIX = new double[][]
{
	{0.0005, 0.0000, 0.0070, 0.0300, 0.1000}, //From flatland to...
	{0.0005, 0.0000, 0.0800, 0.1800, 0.2000}, //From ocean to...
	{0.0001, 0.0000, 0.0050, 0.0500, 0.1000}, //From hills to...
	{0.0001, 0.0000, 0.0010, 0.0100, 0.1000}, //From mountains to...
	{0.0001, 0.0000, 0.0010, 0.0100, 0.1000}  //From peaks to...
};