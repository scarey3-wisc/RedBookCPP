#include "MeshConnection.h"
#include "MeshMidPoint.h"
#include <glm.hpp>
#include <iostream>
#include <unordered_set>
#include "MathToolkit.h"

using namespace std;

MeshConnection::MeshConnection(MeshPoint* first, MeshPoint* second):
	a(first), b(second), mid(nullptr), river(false), ridgeline(false)
{

}

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

vector<MeshPoint*>
MeshConnection::GetQuadCorners()
{

	MeshPoint* start = a;
	MeshPoint* other = b;
	if (b->GetDetailLevel() > a->GetDetailLevel())
	{
		start = b;
		other = a;
	}

	vector<MeshPoint*> quadCorners;
	MeshPoint* startParentA = nullptr;
	MeshPoint* startParentB = nullptr;
	if (MeshMidPoint* startAsMid = dynamic_cast<MeshMidPoint*>(start))
	{
		startParentA = startAsMid->GetParentA();
		startParentB = startAsMid->GetParentB();
	}
	start->ForEachAdjacent([&quadCorners, startParentA, startParentB, other](MeshPoint* candCorner) {
		if (candCorner == startParentA)
			return;
		if (candCorner == startParentB)
			return;
		if (candCorner->AdjacentContains(other))
			quadCorners.push_back(candCorner);
	});
	if (quadCorners.size() > 2)
	{
		cout << "More than 4 corners in a quad is impossible - deploying repair routine!" << endl;
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

MeshMidPoint* 
MeshConnection::GetMidOrTempMid()
{
	if (mid != nullptr)
		return mid;
	if (temporaryMid.has_value())
		temporaryMid->ReassignTemporaryMidPoint(this);
	else
		temporaryMid.emplace(this, false);
	return &*temporaryMid;
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
	std::cout << "";
	MeshConnection* check = a->GetConnection(b);
	if (check == nullptr)
		check = b->GetConnection(a);
	return check;
}


void 
MeshConnection::RemoveIllegalAdjacencies(MeshPoint* s, MeshPoint* t, vector<MeshPoint*>& quadCorners)
{
	glm::dvec2 sG(s->x, s->y);
	glm::dvec2 tG(t->x, t->y);
	unordered_set<MeshPoint*> disqualified;
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

			glm::dvec2 stOne = MathToolkit::GetIntersectionAsST(sG, tG, s_a, t_b);
			glm::dvec2 stTwo = MathToolkit::GetIntersectionAsST(sG, tG, s_b, t_a);
			if (stOne[0] > 0 && stOne[0] < 1 && stOne[1] > 0 && stOne[1] < 1)
			{
				//we found the problem: s->a and t->b cannot simultaneously be connected
				if (s_a.length() < t_b.length())
				{
					//remove the longer t-b edge
					disqualified.insert(b);
					MeshConnection* con = t->GetConnection(b);
					if (con == nullptr)
						con = b->GetConnection(t);
					if (con != nullptr && con->MidInitialized())
					{
						cout << "Trying to remove an invalid edge, but that edge has a midpoint!" << endl;
					}
					t->ForceRemoveAdjacency(b);
					b->ForceRemoveAdjacency(t);
				}
				else
				{
					//remove the longer s-a edge
					disqualified.insert(a);
					MeshConnection* con = s->GetConnection(a);
					if (con == nullptr)
						con = a->GetConnection(s);
					if (con != nullptr && con->MidInitialized())
					{
						cout << "Trying to remove an invalid edge, but that edge has a midpoint!" << endl;
					}
					s->ForceRemoveAdjacency(a);
					a->ForceRemoveAdjacency(s);
				}
			}
			else if (stTwo[0] > 0 && stTwo[0] < 1 && stTwo[1] > 0 && stTwo[1] < 1)
			{
				//we found the problem: s->b and t->a cannot simultaneously be connected
				if (s_b.length() < t_a.length())
				{
					//remove the longer t-a edge
					disqualified.insert(a);
					MeshConnection* con = t->GetConnection(a);
					if (con == nullptr)
						con = a->GetConnection(t);
					if (con != nullptr && con->MidInitialized())
					{
						cout << "Trying to remove an invalid edge, but that edge has a midpoint!" << endl;
					}
					t->ForceRemoveAdjacency(a);
					a->ForceRemoveAdjacency(t);
				}
				else
				{
					//remove the longer s-b edge
					disqualified.insert(b);
					MeshConnection* con = s->GetConnection(b);
					if (con == nullptr)
						con = b->GetConnection(s);
					if (con != nullptr && con->MidInitialized())
					{
						cout << "Trying to remove an invalid edge, but that edge has a midpoint!" << endl;
					}
					s->ForceRemoveAdjacency(b);
					b->ForceRemoveAdjacency(s);
				}
			}
		}
	}
	if (quadCorners.size() - disqualified.size() < 2)
	{
		cout << "Oh boy, we disqualified too many... aborting" << endl;
		return;
	}
	if (quadCorners.size() - disqualified.size() > 2)
	{
		cout << "Oops, we didn't disqualify enough" << endl;
	}
	for (MeshPoint* p : disqualified)
		quadCorners.erase(remove(quadCorners.begin(), quadCorners.end(), p), quadCorners.end());
}