#include "ContinentGenAlgorithms.h"
#include "MeshConnection.h"
#include <iostream>
#include <chrono>

void 
ContinentGenAlgorithms::RunDrainageCalculation(vector<SamplePoint*>& continent)
{
	TectonicUpliftAlgorithm myAlgo;
	myAlgo.InitGraphVertices(continent);
	myAlgo.CalcStreamTrees();
	myAlgo.RunSingleDrainageCalc();
	myAlgo.ResetTemporaryStorage();
}

//------------------------------------------------------------------------------

void 
ContinentGenAlgorithms::BlurUpliftForTectonicAlgorithm(vector<SamplePoint*>& continent, int numIterations, double centerWeight)
{
	vector<double> nextUplift;
	nextUplift.reserve(continent.size());
	for (int i = 0; i < continent.size(); i++)
		nextUplift.push_back(0.);

	for (int i = 0; i < numIterations; i++)
	{
		for (int index = 0; index < continent.size(); index++)
		{
			SamplePoint* sp = continent[index];
			double upliftSum = centerWeight * sp->GetTectonicUplift();
			double weightTotal = centerWeight;
			for (SamplePoint* adj : sp->GetAdjacentSamples())
			{
				if (adj->IsOcean())
					continue;
				weightTotal += 1;
				upliftSum += adj->GetTectonicUplift();
			}
			double newUplift = upliftSum / weightTotal;
			nextUplift[index] = newUplift;
		}

		for (int index = 0; index < continent.size(); index++)
		{
			SamplePoint* sp = continent[index];
			double newUplift = nextUplift[index];
			sp->OverrideTectonicUplift(newUplift);
		}
	}
}

//------------------------------------------------------------------------------

int 
ContinentGenAlgorithms::RunTectonicUpliftAlgorithm(
	vector<SamplePoint*>& continent, double dt, double k,
	int maxIterations, double completionPercentage)
{
	TectonicUpliftAlgorithm myAlgo;
	myAlgo.InitGraphVertices(continent);
	int changeThreshold = (int)(myAlgo.NumIteriorPoints() * completionPercentage);
	int numIterations = 0;
	chrono::steady_clock::time_point currTime, nextTime;
	chrono::steady_clock::duration deltaTime;
	for (; numIterations < maxIterations; numIterations++)
	{
		cout << "It #: " << numIterations << ", ";
		currTime = chrono::steady_clock::now();
		myAlgo.CalcStreamTrees();
		nextTime = chrono::steady_clock::now();
		deltaTime = nextTime - currTime;
		currTime = nextTime;
		cout << "Tree: " << duration_cast<chrono::milliseconds>(deltaTime).count() << ", ";
		myAlgo.UpdateStreamTreesWithLakeOverflow();
		nextTime = chrono::steady_clock::now();
		deltaTime = nextTime - currTime;
		currTime = nextTime;
		cout << "Lake: " << duration_cast<chrono::milliseconds>(deltaTime).count() << ", ";
		myAlgo.RunFullSolveStep(dt, k);
		nextTime = chrono::steady_clock::now();
		deltaTime = nextTime - currTime;
		currTime = nextTime;
		cout << "Solver Total: " << duration_cast<chrono::milliseconds>(deltaTime).count() << ", ";
		int changes = myAlgo.UpdateSinkIndicies();
		nextTime = chrono::steady_clock::now();
		deltaTime = nextTime - currTime;
		currTime = nextTime;
		cout << "Topo: " << duration_cast<chrono::milliseconds>(deltaTime).count() << ", ";
		cout << "Changes: " << changes << "(" << changeThreshold << ")" << endl;
		if (changes < changeThreshold)
			break;
	}

	myAlgo.ResetTemporaryStorage();
	return numIterations;
}

//------------------------------------------------------------------------------

int 
TectonicUpliftAlgorithm::UpdateSinkIndicies()
{
	int count = 0;
	for (int i = 0; i < boundary.size() + interior.size(); i++)
	{
		if (i >= boundary.size())
		{
			size_t chosenIndex = i - boundary.size();
			MeshPoint* mp = interior[chosenIndex];
			int previousSinkIndex = sinkIndices[i];
			int nextSinkIndex = -1;
			if (mp->GetDrainageSink() != nullptr)
				nextSinkIndex = mp->GetDrainageSink()->GetContainerIndex();
			if (previousSinkIndex != nextSinkIndex)
			{
				sinkIndices[i] = nextSinkIndex;
				count++;
			}
		}
		else
		{
			size_t chosenIndex = i;
			MeshPoint* mp = boundary[chosenIndex];
			int previousSinkIndex = sinkIndices[i];
			int nextSinkIndex = -1;
			if (mp->GetDrainageSink() != nullptr)
				nextSinkIndex = mp->GetDrainageSink()->GetContainerIndex();
			if (previousSinkIndex != nextSinkIndex)
			{
				sinkIndices[i] = nextSinkIndex;
				count++;
			}
		}
		
	}
	
	return count;
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::RunSingleDrainageCalc()
{
	int validRoots = 0;

	for (MeshPoint* root : roots)
	{
		if (root->GetDrainageSink() == nullptr)
			continue;
		validRoots++;
		vector<MeshPoint*> tree = GetNodeTreeByDistanceFromRoot(root);
		CalculateDrainageArea(tree);
	}
	cout << "<" + validRoots << ">";
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::RunFullSolveStep(double dt, double k)
{
	int numValidRoots = 0;
	chrono::steady_clock::duration totalBFS = chrono::steady_clock::duration::zero();
	chrono::steady_clock::duration totalDrain = chrono::steady_clock::duration::zero();
	chrono::steady_clock::duration totalSolve = chrono::steady_clock::duration::zero();
	for (MeshPoint* root : roots)
	{
		if (root->GetDrainageSink() != nullptr)
			continue;
		numValidRoots++;

		chrono::steady_clock::time_point start = chrono::steady_clock::now();
		vector<MeshPoint*> tree = GetNodeTreeByDistanceFromRoot(root);
		chrono::steady_clock::time_point treeDone = chrono::steady_clock::now();
		CalculateDrainageArea(tree);
		chrono::steady_clock::time_point drainageDone = chrono::steady_clock::now();
		RunImplicitSolver(tree, dt, k);
		chrono::steady_clock::time_point solveDone = chrono::steady_clock::now();
		totalBFS += (treeDone - start);
		totalDrain += (drainageDone - treeDone);
		totalSolve += (solveDone - drainageDone);
	}

	cout << "(Trees: " << roots.size();
	cout << ", bfs: " << chrono::duration_cast<chrono::milliseconds>(totalBFS).count(); 
	cout << ", drain : " << chrono::duration_cast<chrono::milliseconds>(totalDrain).count();
	cout << ", solve : " << chrono::duration_cast<chrono::milliseconds>(totalSolve).count() << ") ";
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::RunImplicitSolver(vector<MeshPoint*>& orderedNodes, double dt, double k)
{
	for (int i = 0; i < orderedNodes.size(); i++)
	{
		MeshPoint* mp = orderedNodes[i];
		if (mp->GetContainerIndex() < boundary.size())
		{
			if (i != 0)
				cout << "Unexpected result in ContinentGenAlgorithms, TectonicUpliftAlgorithm, RunImplicitSolver!" << endl;
			continue;
		}

		if (mp->IsOcean())
			continue;
		MeshPoint* sink = mp->GetDrainageSink();
		if (sink == nullptr)
		{
			cout << "Very so bad result in ContinentGenAlgorithms, TectonicUpliftAlgorithm, RunImplicitSolver!" << endl;
		}
		
		SamplePoint* sp = dynamic_cast<SamplePoint*>(mp);
		SamplePoint* sinkSp = dynamic_cast<SamplePoint*>(sink);

		if (mp->IsWaterPoint() && sink->IsWaterPoint())
		{
			mp->SetElevation(sink->GetElevation());
			continue;
		}

		double uplift = mp->GetTectonicUplift();

		double dist = mp->DistTo(sink);
		double erosionCoeff = k * pow(mp->GetDrainageArea(), 0.5) / dist;
		double denominator = 1 + erosionCoeff * dt;
		double parenth = uplift + erosionCoeff * sink->GetElevation();
		double numerator = mp->GetElevation() + dt * parenth;
		double result = numerator / denominator;

		double maxResult = mp->GetMaxGrade() * MeshPoint::ConvertVoronoiDistToMeters(dist) + sink->GetElevation();
		if (result > maxResult)
			result = maxResult;
		mp->SetElevation(result);
	}
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::CalculateDrainageArea(vector<MeshPoint*>& orderedNodes)
{
	for (int i = orderedNodes.size() - 1; i >= 0; i--)
	{
		MeshPoint* mp = orderedNodes[i];
		mp->CalculateDrainageArea();
	}
}

//------------------------------------------------------------------------------

vector<MeshPoint*> 
TectonicUpliftAlgorithm::GetNodeTreeByDistanceFromRoot(MeshPoint* root)
{
	vector<MeshPoint*> found;
	deque<MeshPoint*> frontier;
	frontier.push_back(root);
	while (!frontier.empty())
	{
		MeshPoint* mp = frontier.front();
		frontier.pop_front();
		found.push_back(mp);
		for (MeshPoint* source : mp->GetDrainageSources())
			frontier.push_back(source);
	}
	return found;
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::UpdateStreamTreesWithLakeOverflow()
{
	//Apply tags of roots all the way up the stack
	//Is this an opportunity to parallelize each root's stack?
	for (MeshPoint* mp : roots)
	{
		int tag = mp->GetContainerIndex();
		UpdateTreeTags(mp, tag);
	}

	vector<GraphEdge> accepted;
	map<double, deque<GraphEdge>> candEdges;

	//Go through the trees originating from the boundary; find edges out of those trees
	for (MeshPoint* root : roots)
	{
		if (root->GetTag() >= boundary.size())
			continue;
		AddCandidateEdges(root, candEdges);
	}
	while (!candEdges.empty())
	{
		GraphEdge cand = candEdges.begin()->second.front();
		candEdges.begin()->second.pop_front();
		if(candEdges.begin()->second.empty())
			candEdges.erase(candEdges.begin());
		if (!cand.IsCandidateEdge(boundary.size()))
			continue;
		accepted.push_back(cand);
		MeshPoint* lakeRoot = interior[cand.source->GetTag() - boundary.size()];
		int boundaryRootTag = cand.sink->GetTag();
		UpdateTreeTags(lakeRoot, boundaryRootTag);
		AddCandidateEdges(lakeRoot, candEdges);
	}

	for (GraphEdge ge : accepted)
	{
		MeshPoint* lake = ge.source;
		while (lake->GetDrainageSink() != nullptr)
			lake = lake->GetDrainageSink();
		lake->ForceAssignDrainage(ge.sink);
	}
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::AddCandidateEdges(MeshPoint* root, map<double, deque<GraphEdge>>& candEdges)
{
	vector<MeshPoint*> childrenStack;
	childrenStack.push_back(root);
	while (!childrenStack.empty())
	{
		MeshPoint* m = childrenStack.back();
		childrenStack.pop_back();
		for (MeshPoint* src : m->GetDrainageSources())
			childrenStack.push_back(src);

		for (MeshPoint* adj : m->GetDirectlyAdjacent())
		{
			if (adj->GetTag() < 0)
				continue;
			//this adjacency is from a tree rooted in the boundary;
			//adj is not connected to a local minima, so we
			//don't need to work with it
			if (adj->GetTag() < boundary.size())
				continue;

			GraphEdge pass(adj, m);
			candEdges[pass.h].push_front(pass);
		}
	}
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::UpdateTreeTags(MeshPoint* root, int newTag)
{
	vector<MeshPoint*> childrenStack;
	childrenStack.push_back(root);
	while (!childrenStack.empty())
	{
		MeshPoint* m = childrenStack.back();
		childrenStack.pop_back();
		m->SetTag(newTag);
		for (MeshPoint* src : m->GetDrainageSources())
			childrenStack.push_back(src);
	}
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::CalcStreamTrees()
{
	for (MeshPoint* mp : interior)
	{
		mp->ResetDrainage();
		mp->ResetTag();
	}
	for (MeshPoint* mp : boundary)
	{
		mp->ResetDrainage();
		mp->ResetTag();
	}
	for (MeshPoint* mp : interior)
		mp->AssignDrainage();

	roots.clear();
	roots.insert(roots.end(), boundary.begin(), boundary.end());
	for (MeshPoint* mp : interior)
	{
		if(mp->GetDrainageSink() == nullptr)
			roots.push_back(mp);
	}
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::InitGraphVertices(vector<SamplePoint*>& continent)
{
	boundary.clear();
	interior.clear();
	MeshPoint::StartNewSearch();
	vector<vector<MeshPoint*>> detailLevels;
	vector<MeshPoint*> detail0;
	for (SamplePoint* s : continent)
	{
		s->MarkAsReached();
		detail0.push_back(s);
	}
	detailLevels.push_back(std::move(detail0));
	while (true)
	{
		vector<MeshPoint*> nextLevel;
		for (MeshPoint* m : detailLevels.back())
		{
			m->ForEachAdjacent([&nextLevel, &m](MeshPoint* adj) {
				if (!adj->Reached())
					return;
				MeshConnection* con = m->GetConnection(adj);
				if (con == nullptr)
					return;
				MeshPoint* mid = con->GetMid();
				if (mid == nullptr)
					return;
				if (mid->Reached())
					return;
				mid->MarkAsReached();
				nextLevel.push_back(mid);
			});
		}
		if (nextLevel.size() == 0)
			break;
		detailLevels.push_back(std::move(nextLevel));
	}
	for (vector<MeshPoint*> almp : detailLevels)
	{
		for (MeshPoint* p : almp)
		{
			p->InitDirectlyAdjacent();
			p->CalculatePersonalDrainageArea();
			interior.push_back(p);
			for (MeshPoint* adj : p->GetDirectlyAdjacent())
			{
				if (adj->Reached())
					continue;
				adj->MarkAsReached();
				adj->InitDirectlyAdjacent();
				adj->CalculatePersonalDrainageArea();
				boundary.push_back(adj);
			}
		}
	}
	recordedContainerIndices.clear();
	for (size_t i = 0; i < boundary.size(); i++)
	{
		recordedContainerIndices.push_back(boundary[i]->GetContainerIndex());
		boundary[i]->SetContainerIndex(i);
	}
	for (size_t i = 0; i < interior.size(); i++)
	{
		recordedContainerIndices.push_back(interior[i]->GetContainerIndex());
		interior[i]->SetContainerIndex(i + boundary.size());
	}
	sinkIndices.clear();
	sinkIndices.reserve(boundary.size() + interior.size());
	for (int i = 0; i < boundary.size(); i++)
		sinkIndices.push_back(-1);
	for (int i = 0; i < interior.size(); i++)
		sinkIndices.push_back(-1);
}

//------------------------------------------------------------------------------

void 
TectonicUpliftAlgorithm::ResetTemporaryStorage()
{
	for (int i = 0; i < boundary.size(); i++)
	{
		MeshPoint* m = boundary[i];
		m->ResetTag();
		m->ResetDirectlyAdjacent();
		m->SetContainerIndex(recordedContainerIndices[i]);
	}
	for (int i = 0; i < interior.size(); i++)
	{
		MeshPoint* m = interior[i];
		m->ResetTag();
		m->ResetDirectlyAdjacent();
		m->SetContainerIndex(recordedContainerIndices[i + boundary.size()]);
	}
}

//------------------------------------------------------------------------------

int 
TectonicUpliftAlgorithm::NumIteriorPoints()
{
	return (int) interior.size();
}