#include "LocalMap.h"
#include "RegionalMap.h"


LocalMap::LocalMap(int x, int y, RegionalMap* parent) : x(x), y(y), parent(parent), activityFlag(false)
{
	/*String name = "Datamap_" + Integer.toString(x) + "_" + Integer.toString(y);
	String fullName = "DMap(" + Integer.toString(parent.GetWorldX());
	fullName += ";" + String.format("%02d", x) + " , ";
	fullName += Integer.toString(parent.GetWorldY());
	fullName += ";" + String.format("%02d", y) + ")";
	heightmap = new DataImage32Decimal(p.GetDirectory(), RegionalMap.K_HEIGHTMAP_FOLDER_NAME, name, fullName, new HeightQuery(), RedBook.heightmaps);
	watermap = new DataImageByte(p.GetDirectory(), RegionalMap.K_WATERMAP_FOLDER_NAME, name, fullName, new WaterQuery(), RedBook.watermaps);
	rainflowmap = new DataImageInt(p.GetDirectory(), RegionalMap.K_RAINFLOWMAP_FOLDER_NAME, name, fullName, new RainflowQuery(), RedBook.rainflowmaps);
	sedimentmap = new DataImage32Decimal(p.GetDirectory(), RegionalMap.K_SEDIMENTMAP_FOLDER_NAME, name, fullName, new SedimentPerlinQuery(), RedBook.sedimentmaps);
	*/
}






bool 
LocalMap::PrepareForEditing(bool height, bool water, bool rainflow)
{
	return false;
	/*if (height)
	{
		heightmap.RemoveFromManagement();
		heightmap.ForceEditReady(false);
		sedimentmap.RemoveFromManagement();
		sedimentmap.ForceEditReady(false);
	}
	if (water)
	{
		watermap.RemoveFromManagement();
		watermap.ForceEditReady(false);
	}
	if (rainflow)
	{
		rainflowmap.RemoveFromManagement();
		rainflowmap.ForceEditReady(false);
	}
	boolean alreadyEditing = false;
	if (height || water)
		alreadyEditing = RecordWaterHeights(false);
	return alreadyEditing;*/
}
void 
LocalMap::CompleteEditing(bool height, bool water, bool rainflow, bool killWaterHeightsRec)
{
	/*if (height || water)
		FixWaterHeightsFromRecord(killWaterHeightsRec);
	if (height)
	{
		heightmap.SaveAllResolutions(false);
		heightmap.GiveToManagement(RedBook.heightmaps);
		sedimentmap.SaveAllResolutions(false);
		sedimentmap.GiveToManagement(RedBook.sedimentmaps);
	}
	if (water)
	{
		watermap.SaveAllResolutions(false);
		watermap.GiveToManagement(RedBook.watermaps);
	}
	if (rainflow)
	{
		rainflowmap.SaveAllResolutions(false);
		rainflowmap.GiveToManagement(RedBook.rainflowmaps);
	}*/
}
void 
LocalMap::LaplacianErosionIteration(int num)
{
	/*boolean alreadyEditing = PrepareForEditing(true, false, false);
	for (int n = 0; n < num; n++)
	{
		heightmap.InitializePhasedDelta();
		sedimentmap.InitializePhasedDelta();
		for (int i = 0; i <= DataImage.trueDim; i++)
		{
			for (int j = 0; j <= DataImage.trueDim; j++)
			{
				double x = 1.0 * i / DataImage.trueDim;
				double y = 1.0 * j / DataImage.trueDim;

				WatermapValue water = WatermapValue.fromByte(watermap.Get(x, y));
				if (water != WatermapValue.NotWater)
					continue;

				double lap = GetHeightLaplacian(x, y);
				double mpp = 1.0 * METER_DIM / DataImage.trueDim; //meters per pixel

				//accumulate sediment
				if (lap > 0)
				{
					double delta = Switches.LAPLACE_EROSION_DEPOSITION_CONSTANT * lap * mpp * mpp;
					heightmap.PrepareDelta(i, j, delta);
					sedimentmap.PrepareDelta(i, j, delta);
				}
				else if (lap < 0)
				{
					double sedMaxDelta = Switches.LAPLACE_EROSION_SEDIMENT_CONSTANT * lap * mpp * mpp;
					double currentSediment = sedimentmap.Get(i, j);
					if (currentSediment > sedMaxDelta * -1)
					{
						heightmap.PrepareDelta(i, j, sedMaxDelta);
						sedimentmap.PrepareDelta(i, j, sedMaxDelta);
					}
					else
					{
						double percent = -1 * currentSediment / sedMaxDelta;
						double sedimentDelta = currentSediment * -1;
						double bedrockPercent = 1 - percent;
						double bedrockMaxDelta = Switches.LAPLACE_EROSION_ROCK_CONSTANT * lap * mpp * mpp;
						double bedrockDelta = bedrockPercent * bedrockMaxDelta;
						heightmap.PrepareDelta(i, j, sedimentDelta + bedrockDelta);
						sedimentmap.PrepareDelta(i, j, sedimentDelta);
					}
				}
			}
		}
		sedimentmap.FlushDelta();
		heightmap.FlushDelta();
	}
	CompleteEditing(true, false, false, !alreadyEditing);*/
}
void 
LocalMap::SendRandomRainErosion(int numDroplets)
{
	/*HashMap<LocalMap, Boolean> used = new HashMap<LocalMap, Boolean>();
	for (int i = 0; i < numDroplets; i++)
	{
		int perCurr = 100 * i / numDroplets;
		int perNext = 100 * (i + 1) / numDroplets;
		if (perCurr != perNext)
			System.out.println(perCurr + "% through sending the rain");
		WaterDroplet nova = new WaterDroplet(this, 1, used, true);
		boolean okay = true;
		while (okay)
		{
			okay = nova.OneErosionStep();
		}
	}
	for (Entry<LocalMap, Boolean> usedMap : used.entrySet())
	{
		usedMap.getKey().CompleteEditing(true, false, false, !usedMap.getValue());
	}*/
}

double 
LocalMap::GetHeight(int px, int py)
{
	return 0;
	//return heightmap.Get(px, py);
}
//local coordinates - going from 0 to 1 in this local map
double 
LocalMap::GetHeight(double lX, double lY)
{
	return 0;
	//return heightmap.Get(lX, lY);
}
glm::dvec2 
LocalMap::GetHeightGradient(double lX, double lY)
{
	glm::dvec2 grad(0, 0);
	//Vec2 grad = heightmap.GetGradient(lX, lY);
	grad /= METER_DIM;
	return grad;
}
double 
LocalMap::GetHeightLaplacian(double lX, double lY)
{
	double lap = 0;
	//double lap = heightmap.GetLaplacian(lX, lY);
	lap /= (METER_DIM * METER_DIM);
	return lap;
}

int 
LocalMap::GetWorldX()
{
	return parent->GetWorldX() * RegionalMap::DIMENSION + x;
}
int 
LocalMap::GetWorldY()
{
	return parent->GetWorldY() * RegionalMap::DIMENSION + y;
}
LocalMap*
LocalMap::GetNorthWest()
{
	if (GetNorth() != nullptr)
		return GetNorth()->GetWest();
	if (GetWest() != nullptr)
		return GetWest()->GetNorth();
	return nullptr;
}
LocalMap*
LocalMap::GetNorthEast()
{
	if (GetNorth() != nullptr)
		return GetNorth()->GetEast();
	if (GetEast() != nullptr)
		return GetEast()->GetNorth();
	return nullptr;
}
LocalMap* 
LocalMap::GetSouthWest()
{
	if (GetSouth() != nullptr)
		return GetSouth()->GetWest();
	if (GetWest() != nullptr)
		return GetWest()->GetSouth();
	return nullptr;
}
LocalMap*
LocalMap::GetSouthEast()
{
	if (GetSouth() != nullptr)
		return GetSouth()->GetEast();
	if (GetEast() != nullptr)
		return GetEast()->GetSouth();
	return nullptr;
}
LocalMap* 
LocalMap::GetNorth()
{
	if (y - 1 < 0)
	{
		RegionalMap* north = parent->GetNorth();
		if (north == nullptr)
			return nullptr;
		return north->GetLocalMapAt(x, RegionalMap::DIMENSION - 1);
	}
	return parent->GetLocalMapAt(x, y - 1);
}
LocalMap* 
LocalMap::GetEast()
{
	if (x + 1 == RegionalMap::DIMENSION)
	{
		RegionalMap* east = parent->GetEast();
		if (east == nullptr)
			return nullptr;
		return east->GetLocalMapAt(0, y);
	}
	return parent->GetLocalMapAt(x + 1, y);
}
LocalMap*
LocalMap::GetWest()
{
	if (x - 1 < 0)
	{
		RegionalMap* west = parent->GetWest();
		if (west == nullptr)
			return nullptr;
		return west->GetLocalMapAt(RegionalMap::DIMENSION - 1, y);
	}
	return parent->GetLocalMapAt(x - 1, y);
}
LocalMap*
LocalMap::GetSouth()
{
	if (y + 1 == RegionalMap::DIMENSION)
	{
		RegionalMap* south = parent->GetSouth();
		if (south == nullptr)
			return nullptr;
		return south->GetLocalMapAt(x, 0);
	}
	return parent->GetLocalMapAt(x, y + 1);
}