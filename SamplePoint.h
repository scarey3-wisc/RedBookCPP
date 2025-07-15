#pragma once
#include "TerrainType.h"
#include "MeshPoint.h"
#include <glm.hpp>

class RegionalMap;

class SamplePoint : public MeshPoint
{
public:

	SamplePoint(double x, double y, RegionalMap* parent);

	void CalculateAdjacencies();
	std::vector<SamplePoint*> GetAdjacentSamples();

	const std::vector<SamplePoint*>& GetRiverInlets() { return riverInlets; }
	SamplePoint* GetRiverOutlet() { return riverOutlet; }
	bool RiverFlowProcessed() { return riverOutlet != nullptr; }
	bool SetRiverOutlet(SamplePoint* p);
	void SendFlowToOutlet();
	std::vector<SamplePoint*> GrabAdjacentUnassignedRiverFlows(bool sendRiverFlow);
	SamplePoint* GetWayDownhill(bool markOutlet, bool sendRiverFlow);
	void ResetRiver();
	double GetRiverFlow() { return riverFlow; }
	void ForceSetRiverFlow(double flow) { riverFlow = flow; }
	void SendRiverFlow();

	RegionalMap* GetRegionalMap() { return parent; }
	bool NearNorthEdge();
	bool NearSouthEdge();
	bool NearWestEdge();
	bool NearEastEdge();

	bool RidgeAssigned();
	void MakeLake();
	void MakeOcean();
	void AssignVoronoiTerrainType();

	double GetMaxGrade() { return maxGrade; }
	double GetBaseSedimentDepth();
	std::vector<double> GetPerlinElevDiffs();
	uint32_t GetDetailLevel() { return 0; }
	bool IsOcean() { return type.IsTerrainOfType(TerrainTemplate::OCEAN); }
	bool IsInlandLake();
	double GetTectonicUplift();

	void OverrideTectonicUplift(double newUplift) { tectonicUpliftOverride = newUplift; }
	void ResetTectonicUplift() { tectonicUpliftOverride = -1; }


	TerrainType type;
	glm::vec3 myColor;
	bool paintMyColor;

private:

	void Init();

	void AssignMaxGrade(double min, double max);

	void SetPeaks();
	void SetMountains();
	void SetOcean();
	void SetHills();
	void SetLakes();
	void SetFlats();


	RegionalMap* parent;
	double riverFlow;
	SamplePoint* riverOutlet;
	std::vector<SamplePoint*> riverInlets;
	double maxGrade;
	double tectonicUpliftOverride;

public:
	static constexpr int NumPerlinElevDiffs() { return 2; }
};

