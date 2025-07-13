#pragma once
#include <vector>
#include "SamplePoint.h"
#include <array>

class VoronoiAlgorithms
{
public:
	static void EnsureRiversFlowDownhill(std::vector<SamplePoint*>& continent);
	static void IncreaseFractureLevel(std::vector<SamplePoint*>& area);
	static void SetRiverFlow(std::vector<SamplePoint*>& continent);
	static void ConvertCoastalLakeToOcean(std::vector<SamplePoint*>& continentArea);
	static void ConvertSeasToLakes(std::vector<SamplePoint*>& continentCoast, int maxLakeSize);

	static glm::dvec3 BarycentricCoordinatesDY(double x, double y, std::array<MeshPoint*, 3> vertices);
	static glm::dvec3 BarycentricCoordinatesDX(double x, double y, std::array<MeshPoint*, 3> vertices);
	static glm::dvec3 BarycentricCoordinates(double x, double y, std::array<MeshPoint*, 3> vertices);
	static glm::dvec3 BarycentricCoordinates(double x, double y, std::array<SamplePoint*, 3> vertices);

	static double TriangleArea(MeshPoint* a, MeshPoint* b, MeshPoint* c);
	static double TriangleArea(std::array<MeshPoint*, 3> vertices);

	static std::array<MeshPoint*, 3> FindContainingTriangle(double x, double y, SamplePoint* seed);
	static std::array<SamplePoint*, 3> FindContainingSampleTriangle(double x, double y, SamplePoint* seed);

private:
	VoronoiAlgorithms() = delete;

	static double SignedTriangleAreaDY(double x1, double y1, double x2, double y2, double x3, double y3, int whichCoord);
	static double SignedTriangleAreaDX(double x1, double y1, double x2, double y2, double x3, double y3, int whichCoord);
};