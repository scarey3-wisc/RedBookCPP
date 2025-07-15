#pragma once
#include <vector>
#include <array>
#include <glm.hpp>

class SamplePoint;
class MeshPoint;
class TerrainTemplate;

class VoronoiAlgorithms
{
public:
	static void EnsureRiversFlowDownhill(std::vector<SamplePoint*>& continent);
	static void IncreaseFractureLevel(std::vector<SamplePoint*>& area);
	static void SetRiverFlow(std::vector<SamplePoint*>& continent);
	static void ConvertCoastalLakeToOcean(std::vector<SamplePoint*>& continentArea);
	static void ConvertSeasToLakes(std::vector<SamplePoint*>& continentCoast, int maxLakeSize);
	static void ClearHeights(std::vector<SamplePoint*>& coastalPoints);
	static void AssignHeights(std::vector<SamplePoint*>& coastalPoints);
	static void SortByElevation(std::vector<MeshPoint*>& continent);
	static void SortByElevation(std::vector<SamplePoint*>& continent);
	static std::vector<SamplePoint*> FindAllOfType(std::vector<SamplePoint*>& area, TerrainTemplate target);
	static std::vector<std::vector<SamplePoint*>> FindTypeClumps(std::vector<SamplePoint*>& seeds, TerrainTemplate target, int maxNum);
	static std::vector<SamplePoint*> FindAllWithinBoundary(SamplePoint* seed, TerrainTemplate boundaryType);
	static std::vector<SamplePoint*> FindAllWithinBoundary(std::vector<SamplePoint*>& seeds, TerrainTemplate boundaryType);
	static std::vector<SamplePoint*> FindBoundaryPoints(SamplePoint* seed, TerrainTemplate boundaryTo);
	static std::vector<SamplePoint*> FindBoundaryPoints(std::vector<SamplePoint*>& seeds, TerrainTemplate boundaryTo);

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

	static double QueryGrade(TerrainTemplate from, SamplePoint* to);
	static double QueryGrade(TerrainTemplate from, TerrainTemplate to, bool min);
	static int GradeQueryIndex(TerrainTemplate type);

	inline static constexpr double GRADE_MAX_MATRIX[5][5] =
	{
		{0.0040, 0.0000, 0.0250, 0.1800, 0.3000}, //From flatland to...
		{0.0040, 0.0000, 0.1500, 0.5000, 0.5000}, //From ocean to...
		{0.0005, 0.0000, 0.0200, 0.1800, 0.3000}, //From hills to...
		{0.0005, 0.0000, 0.0100, 0.0300, 0.3000}, //From mountains to...
		{0.0005, 0.0000, 0.0100, 0.0100, 0.3000}  //From peaks to...
	};
	inline static constexpr double GRADE_MIN_MATRIX[5][5] =
	{
		{0.0005, 0.0000, 0.0070, 0.0300, 0.1000}, //From flatland to...
		{0.0005, 0.0000, 0.0800, 0.1800, 0.2000}, //From ocean to...
		{0.0001, 0.0000, 0.0050, 0.0500, 0.1000}, //From hills to...
		{0.0001, 0.0000, 0.0010, 0.0100, 0.1000}, //From mountains to...
		{0.0001, 0.0000, 0.0010, 0.0100, 0.1000}  //From peaks to...
	};
};