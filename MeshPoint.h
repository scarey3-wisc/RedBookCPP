#pragma once

#include <vector>
#include <unordered_map>
#include <shared_mutex>
#include <limits>
#include <cmath>

class MeshConnection;

class MeshPoint
{
	friend class MeshConnection;
public:
	MeshPoint(double x, double y) : MeshPoint() 
	{
		this->x = x;
		this->y = y;
	}
	void MarkAsReached() { searchID = CurrentSearchID; }
	bool Reached() { return searchID == CurrentSearchID; }

	void ResetElevation() { currElev = std::numeric_limits<double>::max(); }
	double GetElevation() { return currElev; }
	void SetElevation(double elev) { currElev = elev; }

	double GetDrainageArea() { return totalDrainageArea; }
	double GetPersonalDrainageArea() { return personalDrainageArea; }
	MeshPoint* GetDrainageSink() { return sink; }
	const std::vector<MeshPoint*>& GetDrainageSources() { return sources; }
	void AssignDrainage();
	void CalculateDrainageArea();
	void ForceAssignDrainage(MeshPoint* targetSink);
	void ResetDrainage();
	void CalculatePersonalDrainageArea();

	void ResetDirectlyAdjacent() { directlyAdjacent.clear(); }
	void InitDirectlyAdjacent();
	const std::vector<MeshPoint*>& GetDirectlyAdjacent() { return directlyAdjacent; }
	MeshPoint* GetClosestDirectAdjacency(MeshPoint* adj);

	void ResetAdjacencies();
	MeshConnection* GetConnection(MeshPoint* adj);
	
	template<typename Func>
	void ForEachAdjacent(Func f) 
	{ 
		std::shared_lock lock(adjMutex);
		for (auto p : adjacent)
			f(p.first);
	}

	size_t GetNumAdjacent()
	{
		std::shared_lock lock(adjMutex);
		return adjacent.size();
	}
	bool AdjacentContains(MeshPoint* p);
	//const std::unordered_map<MeshPoint*, MeshConnection*>& GetAdjacent() { return adjacent; }

	bool ContainerIndexSet() { return containerIndex >= 0; }
	void ResetContainerIndex() { containerIndex = -1; }
	int GetContainerIndex() { return containerIndex; }
	void SetContainerIndex(int index) { containerIndex = index; }

	bool TagSet() { return tag >= 0; }
	void ResetTag() { tag = -1; }
	int GetTag() { return tag; }
	void SetTag(int t) { tag = t; }

	double DistTo(MeshPoint* vp) { 
		return sqrt((vp->x - x) * (vp->x - x) + (vp->y - y) * (vp->y - y));
	}

	bool HasZeroElevDiffs();

	virtual double GetMaxGrade() = 0;
	virtual double GetTectonicUplift() = 0;
	virtual bool IsOcean() = 0;
	virtual bool IsInlandLake() = 0;
	bool IsWaterPoint()
	{
		return IsOcean() || IsInlandLake();
	}
	virtual std::vector<double> GetPerlinElevDiffs() = 0;
	virtual uint32_t GetDetailLevel() = 0;


	double x, y;

protected: 
	void ForceRemoveAdjacency(MeshPoint* target);
	void ForceOneWayAdjacency(MeshPoint* target);
	void MarkAdjacent(MeshPoint* p);

private:
	MeshPoint() : x(0), y(0), 
		searchID(0), containerIndex(-1), tag(-1), 
		totalDrainageArea(0), personalDrainageArea(0),
		currElev(0), sink(nullptr){}

	double currElev;
	std::unordered_map<MeshPoint*, MeshConnection*> adjacent;
	std::shared_mutex adjMutex;

	//Temporary storage
	int searchID;
	int containerIndex;
	int tag;

	//for uplift algorithm
	MeshPoint* sink;
	std::vector<MeshPoint*> sources;
	double totalDrainageArea;
	double personalDrainageArea;
	std::vector<MeshPoint*> directlyAdjacent;


public:
	static void StartNewSearch() { ++CurrentSearchID; }
	static double ConvertVoronoiDistToMeters(double d);

private:
	inline static int CurrentSearchID = 0;
};