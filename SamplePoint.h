#pragma once
#include "MeshPoint.h"
#include "TerrainType.h"
#include <glm.hpp>


class SamplePoint : public MeshPoint
{
public:
	TerrainType type;
	glm::vec3 myColor;
	bool paintMyColor;

private:
	RegionalMap parent;
	double riverFlow;
	SamplePoint* riverOutlet;
	std::vector<SamplePoint*> riverInlets;
	double maxGrade;
	double tectonicUpliftOverride;
};

