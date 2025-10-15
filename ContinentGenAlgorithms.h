#pragma once

#include "vector"
#include "map"
#include "deque"
#include "SamplePoint.h"
using namespace std;

struct ContinentGenAlgorithms
{
	static void RunDrainageCalculation(vector<SamplePoint*>& continent);
	static void BlurUpliftForTectonicAlgorithm(vector<SamplePoint*>&, int numIterations, double centerWeight);
	static int RunTectonicUpliftAlgorithm(vector<SamplePoint*>&, double dt, double k, int maxIterations, double completionPercentage);
};

class TectonicUpliftAlgorithm
{
private:
	struct GraphEdge
	{
		MeshPoint* source;
		MeshPoint* sink;
		double h;
		GraphEdge(MeshPoint* source, MeshPoint* sink) : 
			source(source), sink(sink)
		{
			h = max(source->GetElevation(), sink->GetElevation());
		}
		bool IsCandidateEdge(int boundarySize)
		{
			if (source->GetTag() == sink->GetTag())
				return false;
			if (source->GetTag() - boundarySize < 0)
				return false;
			return true;
		}
	};

public:
	int UpdateSinkIndicies();
	void RunSingleDrainageCalc();
	void RunFullSolveStep(double dt, double k);
	void RunImplicitSolver(vector<MeshPoint*>&, double dt, double k);
	void CalculateDrainageArea(vector<MeshPoint*>&);
	vector<MeshPoint*> GetNodeTreeByDistanceFromRoot(MeshPoint* root);
	void UpdateStreamTreesWithLakeOverflow();
	void CalcStreamTrees();
	void InitGraphVertices(vector<SamplePoint*>& continent);
	void ResetTemporaryStorage();
	int NumIteriorPoints();

private:
	void AddCandidateEdges(MeshPoint* root, map<double, deque<GraphEdge>>& candEdges);
	void UpdateTreeTags(MeshPoint* root, int newTag);

private: 
	vector<MeshPoint*> boundary;
	vector <MeshPoint*> interior;
	vector<MeshPoint*> roots;
	vector<int> sinkIndices;
	vector<int> recordedContainerIndices;
};

