#pragma once
#include "MeshMidPoint.h"
#include "MeshPoint.h"
#include <optional>
class MeshConnection
{
public:
	MeshConnection(MeshPoint* first, MeshPoint* second);

	uint8_t GetLargerDetailLevel();
	uint8_t GetSmallerDetailLevel();
	double GetLength() { return a->DistTo(b); }

	bool IsRidgeline() { return ridgeline; }
	bool IsRiver() { return river; }
	void SetRidgeline() { ridgeline = true; }
	void SetRiver() { river = true; }
	void ResetRidgeline() { ridgeline = false; }
	void ResetRiver() { river = false; }

	std::vector<MeshPoint*> GetQuadCorners();

	bool MidInitialized() { return mid != nullptr; }
	MeshMidPoint* GetMid() { return GetMid(false); }
	MeshMidPoint* GetMid(bool spawnIfNeeded);
	MeshMidPoint* GetMidOrTempMid();


	MeshPoint* a;
	MeshPoint* b;

private:
	MeshMidPoint* mid;
	std::optional<MeshMidPoint> temporaryMid;
	bool ridgeline;
	bool river;

public:
	static MeshConnection* FindConnection(MeshPoint* a, MeshPoint* b);

private:
	static void RemoveIllegalAdjacencies(MeshPoint* s, MeshPoint* t, std::vector<MeshPoint*>& quadCorners);
};

