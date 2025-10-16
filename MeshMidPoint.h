#pragma once
#include <vector>
#include "MeshPoint.h"

class MeshConnection;

class MeshMidPoint : public MeshPoint
{
public:
	MeshMidPoint(MeshConnection* parent, bool permanent);
	void ReassignTemporaryMidPoint(MeshConnection* parent);
	MeshPoint* GetParentA() { return a; }
	MeshPoint* GetParentB() { return b; }

	void CopyPerlinElevDiffs(std::vector<double> a) { perlinElevDiffs = a; }
	void AveragePerlinElevDiffs(std::vector<double> a, std::vector<double> b);


	std::vector<double> GetPerlinElevDiffs() { return perlinElevDiffs; }
	uint32_t GetDetailLevel() { return detailLevel; }
	void SetInlandLake() { myWaterType = InlandLake; }
	void SetOcean() { myWaterType = Ocean; }
	void SetNotWater() { myWaterType = NotWater; }
	bool IsOcean() { return myWaterType == Ocean; }
	bool IsInlandLake() { return myWaterType == InlandLake; }
	void SetTectonicUplift(double tu) { tectonicUplift = tu; }
	double GetTectonicUplift() { return tectonicUplift; }
	double GetMaxGrade() { return maxGrade; }

private:
	void InitInterpolation(MeshConnection* parent, bool interpolateElevation);
	void InitConnections(MeshConnection* parent);

	uint8_t myWaterType;
	std::vector<double> perlinElevDiffs;
	uint8_t detailLevel;
	MeshPoint* a;
	MeshPoint* b;
	double tectonicUplift;
	double maxGrade;

	inline static constexpr uint8_t NotWater = 0;
	inline static constexpr uint8_t Ocean = 1;
	inline static constexpr uint8_t InlandLake = 2;
};

