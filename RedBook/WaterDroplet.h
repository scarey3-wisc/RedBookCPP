#pragma once
#include <unordered_map>
class LocalMap;

class WaterDroplet
{
public:
	WaterDroplet(LocalMap* parent, double initWater, std::unordered_map<LocalMap*, bool>& activeList, bool autoLock) :
		parent(parent), water(initWater), activeLocalMaps(activeList) {}
	bool OneErosionStep() { return false; }
private:
	LocalMap* parent;
	double water;
	std::unordered_map<LocalMap*, bool> activeLocalMaps;
};

