#include "MeshMidPoint.h"
#include "MeshConnection.h"

using namespace std;

MeshMidPoint::MeshMidPoint(MeshConnection* parent, bool permanent) : 
	MeshPoint((parent->a->x + parent->b->x) / 2, (parent->a->y + parent->b->y) / 2),
	detailLevel(0), a(nullptr), b(nullptr), tectonicUplift(0), maxGrade(0), myWaterType(0)
{
	a = parent->a;
	b = parent->b;
	myWaterType = NotWater;
	if (!permanent)
	{
		if (a->IsOcean() && b->IsOcean())
			myWaterType = Ocean;
		else if (a->IsInlandLake() && b->IsInlandLake())
			myWaterType = InlandLake;
		else if (a->IsOcean() && b->IsInlandLake())
			myWaterType = InlandLake;
		else if (b->IsInlandLake() && a->IsOcean())
			myWaterType = InlandLake;
	}
	InitInterpolation(parent, true);

	if (!permanent)
		return;
	InitConnections(parent);
}
/*public MeshMidPoint(MeshConnection parent, DataInputStream dis)
{
	super(dis);
	a = parent.a;
	b = parent.b;
	InitInterpolation(parent, false);
	try
	{
		myWaterType = WaterType.fromByte(dis.readByte());
		for (int i = 0; i < perlinElevDiffs.length; i++)
		{
			perlinElevDiffs[i] = dis.readDouble();
		}
	}
	catch (IOException e)
	{
		e.printStackTrace();
	}

	InitConnections(parent);
}
public MeshMidPoint(MeshConnection parent, Iterator<String> tokenStream)
{
	super(tokenStream);
	a = parent.a;
	b = parent.b;
	InitInterpolation(parent, false);
	myWaterType = WaterType.fromByte(Byte.parseByte(tokenStream.next()));
	for (int i = 0; i < perlinElevDiffs.length; i++)
	{
		perlinElevDiffs[i] = Double.parseDouble(tokenStream.next());
	}
	InitConnections(parent);
}
public static boolean ConsumeDescription(DataInputStream dis)
{
	try {
		MeshPoint.ConsumeDescription(dis);
		dis.readByte();
		for (int i = 0; i < SamplePoint.NumPerlinElevDiffs(); i++)
		{
			dis.readDouble();
		}
		return true;
	}
	catch (IOException e) {
		return false;
	}
}*/

void 
MeshMidPoint::InitInterpolation(MeshConnection* parent, bool interpolateElevation)
{
	detailLevel = (uint8_t)(parent->GetLargerDetailLevel() + 1);
	AveragePerlinElevDiffs(parent->a->GetPerlinElevDiffs(), parent->b->GetPerlinElevDiffs());
	if (interpolateElevation)
		SetElevation(a->GetElevation() / 2 + b->GetElevation() / 2);
	tectonicUplift = a->GetTectonicUplift() / 2 + b->GetTectonicUplift() / 2;
	maxGrade = a->GetMaxGrade() / 2 + b->GetMaxGrade() / 2;
}
void 
MeshMidPoint::InitConnections(MeshConnection* parent)
{
	ForceOneWayAdjacency(a);
	ForceOneWayAdjacency(b);
	MeshConnection* connA = GetConnection(a);
	MeshConnection* connB = GetConnection(b);
	if (parent->IsRidgeline())
	{
		connA->SetRidgeline();
		connB->SetRidgeline();
	}
	if (parent->IsRiver())
	{
		connA->SetRiver();
		connB->SetRiver();
	}
	for (MeshPoint* c : parent->GetQuadCorners())
	{
		MeshConnection* connOne = MeshConnection::FindConnection(a, c);
		if (connOne != nullptr && connOne->MidInitialized())
			MarkAdjacent(connOne->GetMid());

		MeshConnection* connTwo = MeshConnection::FindConnection(b, c);
		if (connTwo != nullptr && connTwo->MidInitialized())
			MarkAdjacent(connTwo->GetMid());
	}
}

/*public String GetDescription()
{
	String desc = super.GetDescription() + " " + Byte.toString(myWaterType.toByte()) + " ";
	for (int i = 0; i < perlinElevDiffs.length; i++)
	{
		double d = perlinElevDiffs[i];
		desc += Double.toString(d);
		if (i != perlinElevDiffs.length - 1)
			desc += " ";
	}
	return desc;
}

public boolean WriteDescription(DataOutputStream dos)
{
	if (!super.WriteDescription(dos))
		return false;
	try {
		dos.writeByte(myWaterType.toByte());
		for (double d : perlinElevDiffs)
			dos.writeDouble(d);
		return true;
	}
	catch (IOException e) {
		return false;
	}
}*/

void 
MeshMidPoint::AveragePerlinElevDiffs(vector<double> a, vector<double> b)
{
	perlinElevDiffs.clear();
	for (int i = 0; i < a.size() && i < b.size(); i++)
		perlinElevDiffs.push_back((a[i] + b[i]) / 2);
}
