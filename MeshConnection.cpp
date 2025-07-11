#include "MeshConnection.h"
#include <glm.hpp>

uint8_t
MeshConnection::GetLargerDetailLevel()
{
	if (a->GetDetailLevel() > b->GetDetailLevel())
		return a->GetDetailLevel();
	return b->GetDetailLevel();
}



uint8_t
MeshConnection::GetSmallerDetailLevel()
{
	if (a->GetDetailLevel() < b->GetDetailLevel())
		return a->GetDetailLevel();
	return b->GetDetailLevel();

}

std::vector<MeshPoint*>
MeshConnection::GetQuadCorners()
{

	MeshPoint* start = a;
	MeshPoint* other = b;
	if (b->GetDetailLevel() > a->GetDetailLevel())
	{
		start = b;
		other = a;
	}

	std::vector<MeshPoint*> quadCorners;
	MeshPoint* startParentA = nullptr;
	MeshPoint* startParentB = nullptr;
	if (start instanceof MeshMidPoint)
	{
		MeshMidPoint* startAsMid = (MeshMidPoint*) start;
		startParentA = startAsMid.GetParentA();
		startParentB = startAsMid.GetParentB();
	}
	for (MeshPoint* candCorner : start.GetAdjacent())
	{
		if (candCorner == startParentA)
			continue;
		if (candCorner == startParentB)
			continue;
		if (candCorner.GetAdjacent().contains(other))
			quadCorners.add(candCorner);
	}
	if (quadCorners.size() > 2)
	{
		System.out.println("More than 4 corners in a quad is impossible - deploying repair routine!");
		RemoveIllegalAdjacencies(start, other, quadCorners);
	}

	return quadCorners;
}

MeshMidPoint* 
MeshConnection::GetMid(bool spawnIfNeeded)
{
	if (mid == nullptr && spawnIfNeeded)
		mid = new MeshMidPoint(this, true);
	return mid;
}
/*public MeshMidPoint GetMid(DataInputStream dis)
{
	if (mid != null)
		MeshMidPoint.ConsumeDescription(dis);
	else
		mid = new MeshMidPoint(this, dis);
	return mid;
}
public MeshMidPoint GetMid(Iterator<String> tokenStream)
{
	MeshMidPoint nova = new MeshMidPoint(this, tokenStream);
	if (mid == null)
		mid = nova;
	return mid;
}*/


MeshConnection* 
MeshConnection::FindConnection(MeshPoint* a, MeshPoint* b)
{
	MeshConnection* check = a->GetConnection(b);
	if (check == nullptr)
		check = b->GetConnection(a);
	return check;
}


void 
MeshConnection::RemoveIllegalAdjacencies(MeshPoint* s, MeshPoint* t, std::vector<MeshPoint*>& quadCorners)
{
	glm::dvec2 sG(s->x, s->y);
	glm::dvec2 tG(t->x, t->y);
	HashSet<MeshPoint> disqualified = new HashSet<MeshPoint>();
	for (int i = 0; i < quadCorners.size(); i++)
	{
		for (int j = i + 1; j < quadCorners.size(); j++)
		{
			MeshPoint* a = quadCorners[i];
			MeshPoint* b = quadCorners[j];
			glm::dvec2 s_a(a->x - s->x, a->y - s->y);
			glm::dvec2 s_b(b->x - s->x, b->y - s->y);
			glm::dvec2 t_a(a->x - t->x, a->y - t->y);
			glm::dvec2 t_b(b->x - t->x, b->y - t->y);

			double[] stOne = Vec2.GetIntersectionAsST(sG, tG, s_a, t_b);
			double[] stTwo = Vec2.GetIntersectionAsST(sG, tG, s_b, t_a);
			if (stOne != null && stOne[0] > 0 && stOne[0] < 1 && stOne[1] > 0 && stOne[1] < 1)
			{
				//we found the problem: s->a and t->b cannot simultaneously be connected
				if (s_a.Len() < t_b.Len())
				{
					//remove the longer t-b edge
					disqualified.add(b);
					MeshConnection con = t.GetConnection(b);
					if (con == null)
						con = b.GetConnection(t);
					if (con != null && con.MidInitialized())
					{
						System.out.println("Trying to remove an invalid edge, but that edge has a midpoint!");
					}
					t.ForceRemoveAdjacency(b);
					b.ForceRemoveAdjacency(t);
				}
				else
				{
					//remove the longer s-a edge
					disqualified.add(a);
					MeshConnection con = s.GetConnection(a);
					if (con == null)
						con = a.GetConnection(s);
					if (con != null && con.MidInitialized())
					{
						System.out.println("Trying to remove an invalid edge, but that edge has a midpoint!");
					}
					s.ForceRemoveAdjacency(a);
					a.ForceRemoveAdjacency(s);
				}
			}
			else if (stTwo != null && stTwo[0] > 0 && stTwo[0] < 1 && stTwo[1] > 0 && stTwo[1] < 1)
			{
				//we found the problem: s->b and t->a cannot simultaneously be connected
				if (s_b.Len() < t_a.Len())
				{
					//remove the longer t-a edge
					disqualified.add(a);
					MeshConnection con = t.GetConnection(a);
					if (con == null)
						con = a.GetConnection(t);
					if (con != null && con.MidInitialized())
					{
						System.out.println("Trying to remove an invalid edge, but that edge has a midpoint!");
					}
					t.ForceRemoveAdjacency(a);
					a.ForceRemoveAdjacency(t);
				}
				else
				{
					//remove the longer s-b edge
					disqualified.add(b);
					MeshConnection con = s.GetConnection(b);
					if (con == null)
						con = b.GetConnection(s);
					if (con != null && con.MidInitialized())
					{
						System.out.println("Trying to remove an invalid edge, but that edge has a midpoint!");
					}
					s.ForceRemoveAdjacency(b);
					b.ForceRemoveAdjacency(s);
				}
			}
		}
	}
	if (quadCorners.size() - disqualified.size() < 2)
	{
		System.out.println("Oh boy, we disqualified too many... aborting");
		return;
	}
	if (quadCorners.size() - disqualified.size() > 2)
	{
		System.out.println("Oops, we didn't disqualify enough");
	}
	for (MeshPoint p : disqualified)
		quadCorners.remove(p);

}