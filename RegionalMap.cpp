#include "RegionalMap.h"
#include "Perlin.h"
#include "SamplePoint.h"
#include "VoronoiAlgorithms.h"
#include "MeshConnection.h"
#include "Switches.h"
#include "WaterDroplet.h"
#include "WorldMap.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

using namespace std;

RegionalMap::RegionalMap(int x, int y, WorldMap* p) : x(x), y(y), parent(p), readyToRender(false)
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			topography[i * DIMENSION + j] = new LocalMap(i, j, this);
		}
	}
	//EnsureBasicDirectoryStructureExists();
	//SaveBasicDescription();
}
/*public RegionalMap(String directory, WorldMap parent)
{
	this.parent = parent;
	String dir = parent.GetDirectory() + File.separator;
	dir += WorldMap.K_REGIONS_FOLDER_NAME + File.separator;
	dir += directory;
	String name = dir + File.separator + "BasicDesc.txt";
	File basicData = new File(name);
	try {
		Scanner std = new Scanner(basicData);
		this.x = std.nextInt();
		this.y = std.nextInt();
		std.nextLine();
		std.close();
	}
	catch (IOException e) {
		return;
	}

	topography = new LocalMap[DIMENSION * DIMENSION];
	for (int i = 0; i < DIMENSION; i++)
		for (int j = 0; j < DIMENSION; j++)
		{
			topography[i * DIMENSION + j] = new LocalMap(i, j, this);
		}
	terrainCells = new SamplePoint[VORONOI_DIM * VORONOI_DIM];
	voronoiList = new ArrayList<SamplePoint>();
	EnsureBasicDirectoryStructureExists();
	readyToRender = false;
}
private String GetRiverFileName(boolean binary)
{
	if (binary)
		return GetDirectory() + File.separator + "RiverDetails.bin";
	else
		return GetDirectory() + File.separator + "RiverDetails.txt";
}
public boolean SaveRiverFile(boolean binary)
{
	File riverFile = new File(GetRiverFileName(binary));
	if (riverFile.exists())
	{
		if (!riverFile.delete())
			return false;
	}
	if (binary)
		return SaveRiverFileBinary(riverFile);
	else
		return SaveRiverFileVerbose(riverFile);
}
private boolean SaveRiverFileBinary(File file)
{
	try {
		DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		dos.writeInt(voronoiList.size());
		for (SamplePoint sp : voronoiList)
		{
			dos.writeDouble(sp.GetRiverFlow());
			if (sp.RiverFlowProcessed())
			{
				dos.writeByte(1);
				sp.GetRiverOutlet().WriteGlobalIdentifier(dos);
			}
			else
				dos.writeByte(0);
		}
		dos.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean SaveRiverFileVerbose(File file)
{
	try
	{
		BufferedWriter wr = new BufferedWriter(new FileWriter(file));
		for (SamplePoint sp : voronoiList)
		{
			String desc = Double.toString(sp.GetRiverFlow());
			desc += " ";
			if (sp.RiverFlowProcessed())
			{
				desc += "1" + " ";
				desc += sp.GetRiverOutlet().GetGlobalIdentifier();
			}
			else
			{
				desc += "0";
			}
			wr.write(desc);
			wr.newLine();
		}
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public boolean LoadRiverFile()
{
	boolean binary = true;
	File riverFile = new File(GetRiverFileName(binary));
	if (!riverFile.exists() || !riverFile.isFile())
		binary = false;
	riverFile = new File(GetRiverFileName(binary));
	if (!riverFile.exists() || !riverFile.isFile())
		return false;
	if (binary)
		return LoadRiverFileBinary(riverFile);
	else
		return LoadRiverFileVerbose(riverFile);
}
private boolean LoadRiverFileBinary(File file)
{
	try {
		DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		int num = dis.readInt();
		boolean hadIssue = false;
		for (int count = 0; count < num; count++)
		{
			SamplePoint p = voronoiList.get(count);
			double flow = dis.readDouble();
			p.ForceSetRiverFlow(flow);
			byte flag = dis.readByte();
			if (flag == 1)
			{
				SamplePoint outlet = SamplePoint.LoadFromGlobalIdentifier(dis, parent);
				p.SetRiverOutlet(outlet);
			}
		}
		dis.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadRiverFileVerbose(File file)
{
	try {
		boolean hadIssue = false;
		Scanner std = new Scanner(file);
		int index = 0;
		while (std.hasNextLine())
		{
			SamplePoint curr = voronoiList.get(index);
			String desc = std.nextLine();
			String[] dets = desc.split(" ");
			double flow = Double.parseDouble(dets[0]);
			curr.ForceSetRiverFlow(flow);
			if (dets[1].equals("1"))
			{
				SamplePoint outlet = SamplePoint.FindFromGlobalIdentifier(dets[2], parent);
				curr.SetRiverOutlet(outlet);
			}
			index++;
		}
		std.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private String GetMidpointFileName(boolean binary)
{
	if (binary)
		return GetDirectory() + File.separator + "MidpointTree.bin";
	else
		return GetDirectory() + File.separator + "MidpointTree.txt";
}
public boolean SaveMidpointTree(boolean binary)
{
	File treeFile = new File(GetMidpointFileName(binary));
	if (treeFile.exists())
	{
		if (!treeFile.delete())
			return false;
	}
	if (binary)
		return SaveMidpointTreeBinary(treeFile);
	else
		return SaveMidpointTreeVerbose(treeFile);
}
private boolean SaveMidpointTreeBinary(File file)
{
	try
	{
		DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		MeshPoint.StartNewSearch();
		ArrayList<ArrayList<MeshMidPoint>> detailLevels = new ArrayList<ArrayList<MeshMidPoint>>();
		ArrayList< ? extends MeshPoint> currentParents = voronoiList;
		for (MeshPoint p : currentParents)
			p.MarkAsReached();
		while (true)
		{
			ArrayList<MeshMidPoint> found = new ArrayList<MeshMidPoint>();
			for (MeshPoint par : currentParents)
			{
				for (MeshPoint adj : par.GetAdjacent())
				{
					MeshConnection con = par.GetConnection(adj);
					if (con == null)
						continue;
					if (!con.MidInitialized())
						continue;
					MeshMidPoint mid = con.GetMid();
					if (mid.Reached())
						continue;
					mid.MarkAsReached();
					mid.SetContainerIndex(found.size());
					found.add(mid);
				}
			}
			if (found.size() == 0)
				break;
			currentParents = found;
			detailLevels.add(found);
		}
		boolean hadProblem = false;
		for (ArrayList<MeshMidPoint> level : detailLevels)
		{
			dos.writeInt(level.size());
			for (MeshMidPoint p : level)
			{
				if (p.GetParentA() instanceof SamplePoint)
				{
					SamplePoint a = (SamplePoint)p.GetParentA();
					dos.writeByte(0);
					hadProblem = hadProblem || !a.WriteGlobalIdentifier(dos);
				}
				else
				{
					dos.writeByte(p.GetParentA().GetDetailLevel());
					dos.writeInt(p.GetParentA().GetContainerIndex());
				}
				if (p.GetParentB() instanceof SamplePoint)
				{
					SamplePoint b = (SamplePoint)p.GetParentB();
					dos.writeByte(0);
					hadProblem = hadProblem || !b.WriteGlobalIdentifier(dos);
				}
				else
				{
					dos.writeByte(p.GetParentB().GetDetailLevel());
					dos.writeInt(p.GetParentB().GetContainerIndex());
				}
				hadProblem = hadProblem || !p.WriteDescription(dos);
			}
		}
		for (ArrayList<MeshMidPoint> level : detailLevels)
		{
			for (MeshMidPoint p : level)
				p.ResetContainerIndex();
		}
		dos.writeInt(0);
		dos.close();
		return !hadProblem;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean SaveMidpointTreeVerbose(File file)
{
	try
	{
		BufferedWriter wr = new BufferedWriter(new FileWriter(file));
		MeshPoint.StartNewSearch();
		ArrayList<ArrayList<MeshMidPoint>> detailLevels = new ArrayList<ArrayList<MeshMidPoint>>();
		ArrayList< ? extends MeshPoint> currentParents = voronoiList;
		for (MeshPoint p : currentParents)
			p.MarkAsReached();
		while (true)
		{
			ArrayList<MeshMidPoint> found = new ArrayList<MeshMidPoint>();
			for (MeshPoint par : currentParents)
			{
				for (MeshPoint adj : par.GetAdjacent())
				{
					MeshConnection con = par.GetConnection(adj);
					if (con == null)
						continue;
					if (!con.MidInitialized())
						continue;
					MeshMidPoint mid = con.GetMid();
					if (mid.Reached())
						continue;
					mid.MarkAsReached();
					mid.SetContainerIndex(found.size());
					found.add(mid);
				}
			}
			if (found.size() == 0)
				break;
			currentParents = found;
			detailLevels.add(found);
		}
		for (ArrayList<MeshMidPoint> level : detailLevels)
		{
			wr.write(Integer.toString(level.size()));
			wr.newLine();
			for (MeshMidPoint p : level)
			{
				String desc = "";
				if (p.GetParentA() instanceof SamplePoint)
				{
					SamplePoint a = (SamplePoint)p.GetParentA();
					desc += "0" + " " + a.GetGlobalIdentifier();
				}
				else
				{
					desc += Byte.toString(p.GetParentA().GetDetailLevel());
					desc += " " + p.GetParentA().GetContainerIndex();
				}
				desc += " ";
				if (p.GetParentB() instanceof SamplePoint)
				{
					SamplePoint b = (SamplePoint)p.GetParentB();
					desc += "0" + " " + b.GetGlobalIdentifier();
				}
				else
				{
					desc += Byte.toString(p.GetParentB().GetDetailLevel());
					desc += " " + p.GetParentB().GetContainerIndex();
				}
				desc += " ";
				desc += p.GetDescription();
				wr.write(desc);
				wr.newLine();
			}
		}
		for (ArrayList<MeshMidPoint> level : detailLevels)
		{
			for (MeshMidPoint p : level)
				p.ResetContainerIndex();
		}
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public boolean LoadMidpointTree()
{
	boolean binary = true;
	File treeFile = new File(GetMidpointFileName(binary));
	if (!treeFile.exists() || !treeFile.isFile())
		binary = false;
	treeFile = new File(GetMidpointFileName(binary));
	if (!treeFile.exists() || !treeFile.isFile())
		return false;
	if (binary)
		return LoadMidpointTreeBinary(treeFile);
	else
		return LoadMidpointTreeVerbose(treeFile);
}
private boolean LoadMidpointTreeBinary(File file)
{
	try {
		boolean hadIssue = false;
		DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		ArrayList<ArrayList<MeshMidPoint>> detailLevels = new ArrayList<ArrayList<MeshMidPoint>>();
		while (true)
		{
			int levelSize = dis.readInt();
			if (levelSize == 0)
				break;
			ArrayList<MeshMidPoint> myLevel = new ArrayList<MeshMidPoint>(levelSize);
			for (int i = 0; i < levelSize; i++)
			{
				int detA = dis.readByte();
				MeshPoint a = null;
				if (detA == 0)
					a = SamplePoint.LoadFromGlobalIdentifier(dis, parent);
				else
				{
					int index = dis.readInt();
					if (index != -1)
						a = detailLevels.get(detA - 1).get(index);
				}
				int detB = dis.readByte();
				MeshPoint b = null;
				if (detB == 0)
					b = SamplePoint.LoadFromGlobalIdentifier(dis, parent);
				else
				{
					int index = dis.readInt();
					if (index != -1)
						b = detailLevels.get(detB - 1).get(index);
				}
				if (a == null || b == null)
				{
					MeshMidPoint.ConsumeDescription(dis);
					myLevel.add(null);
					continue;
				}
				MeshConnection con = a.GetConnection(b);
				if (con == null)
					con = b.GetConnection(a);
				if (con == null)
				{
					MeshMidPoint.ConsumeDescription(dis);
					myLevel.add(null);
					continue;
				}
				MeshMidPoint mid = con.GetMid(dis);
				myLevel.add(mid);
			}
			detailLevels.add(myLevel);
		}
		dis.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadMidpointTreeVerbose(File file)
{
	try {
		boolean hadIssue = false;
		Scanner std = new Scanner(file);
		ArrayList<ArrayList<MeshMidPoint>> detailLevels = new ArrayList<ArrayList<MeshMidPoint>>();
		while (true)
		{
			if (!std.hasNextInt())
				break;
			int levelSize = std.nextInt();
			std.nextLine();
			ArrayList<MeshMidPoint> myLevel = new ArrayList<MeshMidPoint>(levelSize);
			for (int i = 0; i < levelSize; i++)
			{
				String desc = std.nextLine();
				String[] dets = desc.split(" ");
				int detA = Integer.parseInt(dets[0]);
				MeshPoint a = null;
				if (detA == 0)
					a = SamplePoint.FindFromGlobalIdentifier(dets[1], parent);
				else
				{
					int index = Integer.parseInt(dets[1]);
					if (index != -1)
						a = detailLevels.get(detA - 1).get(index);
				}
				int detB = Integer.parseInt(dets[2]);
				MeshPoint b = null;
				if (detB == 0)
					b = SamplePoint.FindFromGlobalIdentifier(dets[3], parent);
				else
				{
					int index = Integer.parseInt(dets[3]);
					if (index != -1)
						b = detailLevels.get(detB - 1).get(index);
				}
				if (a == null || b == null)
				{
					myLevel.add(null);
					continue;
				}
				MeshConnection con = a.GetConnection(b);
				if (con == null)
					con = b.GetConnection(a);
				if (con == null)
				{
					myLevel.add(null);
					continue;
				}
				Iterator<String> tokenStream = Arrays.stream(Arrays.copyOfRange(dets, 4, dets.length)).iterator();
				MeshMidPoint mid = con.GetMid(tokenStream);
				myLevel.add(mid);
			}
			detailLevels.add(myLevel);
		}
		std.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private String GetSampleAdjFileName(boolean binary)
{
	if (binary)
		return GetDirectory() + File.separator + "SampleAdjacencies.bin";
	else
		return GetDirectory() + File.separator + "SampleAdjacencies.txt";
}
public boolean SaveSampleAdjacencies(boolean binary)
{
	File sampleAdjFile = new File(GetSampleAdjFileName(binary));
	if (sampleAdjFile.exists())
	{
		if (!sampleAdjFile.delete())
			return false;
	}
	if (binary)
		return SaveSampleAdjacenciesBinary(sampleAdjFile);
	else
		return SaveSampleAdjacenciesVerbose(sampleAdjFile);
}
private boolean SaveSampleAdjacenciesBinary(File file)
{
	try {
		DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		dos.writeInt(voronoiList.size());
		for (SamplePoint sp : voronoiList)
			sp.WriteAdjacencies(dos);
		dos.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean SaveSampleAdjacenciesVerbose(File file)
{
	try
	{
		BufferedWriter wr = new BufferedWriter(new FileWriter(file));
		for (SamplePoint sp : voronoiList)
		{
			String desc = sp.GetAdjacencyDescription();
			wr.write(desc);
			wr.newLine();
		}
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public boolean LoadSampleAdjacencies()
{
	boolean binary = true;
	File sampleAdjFile = new File(GetSampleAdjFileName(binary));
	if (!sampleAdjFile.exists() || !sampleAdjFile.isFile())
		binary = false;
	sampleAdjFile = new File(GetSampleAdjFileName(binary));
	if (!sampleAdjFile.exists() || !sampleAdjFile.isFile())
		return false;
	if (binary)
		return LoadSampleAdjacenciesBinary(sampleAdjFile);
	else
		return LoadSampleAdjacenciesVerbose(sampleAdjFile);
}
private boolean LoadSampleAdjacenciesBinary(File file)
{
	try {
		DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		int num = dis.readInt();
		boolean hadIssue = false;
		for (int count = 0; count < num; count++)
		{
			SamplePoint p = voronoiList.get(count);
			if (!p.LoadAdjacencies(dis, parent))
				hadIssue = true;
		}
		dis.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadSampleAdjacenciesVerbose(File file)
{
	try {
		boolean hadIssue = false;
		Scanner std = new Scanner(file);
		int index = 0;
		while (std.hasNextLine())
		{
			String desc = std.nextLine();
			SamplePoint curr = voronoiList.get(index);
			curr.SetAdjacenciesFromDescription(desc, parent);
			index++;
		}
		std.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private String GetSampleListFileName(boolean binary)
{
	if (binary)
		return GetDirectory() + File.separator + "SampleList.bin";
	else
		return GetDirectory() + File.separator + "SampleList.txt";
}
public boolean SaveSampleList(boolean binary)
{
	File sampleListFile = new File(GetSampleListFileName(binary));
	if (sampleListFile.exists())
	{
		if (!sampleListFile.delete())
			return false;
	}
	if (binary)
		return SaveSampleListBinary(sampleListFile);
	else
		return SaveSampleListVerbose(sampleListFile);

}
private boolean SaveSampleListBinary(File file)
{
	try {
		DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
		dos.writeInt(voronoiList.size());
		for (SamplePoint sp : voronoiList)
			sp.WriteDescription(dos);
		dos.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean SaveSampleListVerbose(File file)
{
	try
	{
		BufferedWriter wr = new BufferedWriter(new FileWriter(file));
		for (SamplePoint sp : voronoiList)
		{
			String desc = sp.GetDescription();
			wr.write(desc);
			wr.newLine();
		}
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
public boolean LoadSampleList()
{
	boolean binary = true;
	File sampleListFile = new File(GetSampleListFileName(binary));
	if (!sampleListFile.exists() || !sampleListFile.isFile())
		binary = false;
	sampleListFile = new File(GetSampleListFileName(binary));
	if (!sampleListFile.exists() || !sampleListFile.isFile())
		return false;
	if (binary)
		return LoadSampleListBinary(sampleListFile);
	else
		return LoadSampleListVerbose(sampleListFile);

}
private boolean LoadSampleListBinary(File file)
{
	if (voronoiList.size() > 0)
		return false;
	try {
		DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		int num = dis.readInt();
		boolean hadIssue = false;
		voronoiList.ensureCapacity(num);
		for (int count = 0; count < num; count++)
		{
			SamplePoint p = new SamplePoint(dis, this);
			if (p.x - GetWorldX() < 0 || p.y - GetWorldY() < 0)
			{
				hadIssue = true;
				continue;
			}
			int i = (int)((p.x - GetWorldX()) * VORONOI_DIM);
			int j = (int)((p.y - GetWorldY()) * VORONOI_DIM);
			if (i >= VORONOI_DIM || j >= VORONOI_DIM)
			{
				hadIssue = true;
				continue;
			}
			if (terrainCells[i * VORONOI_DIM + j] == null)
			{
				terrainCells[i * VORONOI_DIM + j] = p;
				int index = voronoiList.size();
				voronoiList.add(p);
				p.SetContainerIndex(index);
			}
		}
		dis.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean LoadSampleListVerbose(File file)
{
	if (voronoiList.size() > 0)
		return false;
	try {
		boolean hadIssue = false;
		Scanner std = new Scanner(file);
		while (std.hasNextLine())
		{
			String desc = std.nextLine();
			Iterator<String> tokenStream = Arrays.stream(desc.split(" ")).iterator();
			SamplePoint p = new SamplePoint(tokenStream, this);
			if (p.x - GetWorldX() < 0 || p.y - GetWorldY() < 0)
			{
				hadIssue = true;
				continue;
			}
			int i = (int)((p.x - GetWorldX()) * VORONOI_DIM);
			int j = (int)((p.y - GetWorldY()) * VORONOI_DIM);
			if (i >= VORONOI_DIM || j >= VORONOI_DIM)
			{
				hadIssue = true;
				continue;
			}
			if (terrainCells[i * VORONOI_DIM + j] == null)
			{
				terrainCells[i * VORONOI_DIM + j] = p;
				int index = voronoiList.size();
				voronoiList.add(p);
				p.SetContainerIndex(index);
			}
		}
		std.close();
		return !hadIssue;
	}
	catch (IOException e) {
		return false;
	}
}
private String GetBasicDescName()
{
	return GetDirectory() + File.separator + "BasicDesc.txt";
}
private boolean SaveBasicDescription()
{
	File basicData = new File(GetBasicDescName());
	if (basicData.exists())
	{
		if (!basicData.delete())
			return false;
	}
	try {
		BufferedWriter wr = new BufferedWriter(new FileWriter(basicData));
		wr.write(Integer.toString(x));
		wr.write(" ");
		wr.write(Integer.toString(y));
		wr.newLine();
		wr.close();
		return true;
	}
	catch (IOException e) {
		return false;
	}
}
private boolean EnsureDirectoryStructureForDataImages(String imageDirectoryPrefix)
{
	for (int dim : DataImage.dimRange)
	{
		String folderName = GetDirectory();
		folderName += File.separator;
		folderName += imageDirectoryPrefix;
		folderName += Integer.toString(dim);
		File lmDir = new File(folderName);
		if (lmDir.exists() && !lmDir.isDirectory())
			lmDir.delete();
		if (!lmDir.exists())
			lmDir.mkdir();
		if (!lmDir.exists() || !lmDir.isDirectory())
			return false;
		if (Switches.CLEAR_IMAGE_CACHES)
		{
			for (File exist : lmDir.listFiles())
			{
				exist.delete();
			}
		}
	}
	return true;
}
private boolean EnsureBasicDirectoryStructureExists()
{
	File topDir = new File(GetDirectory());
	if (topDir.exists() && !topDir.isDirectory())
		topDir.delete();
	if (!topDir.exists())
		topDir.mkdirs();
	if (!topDir.exists() || !topDir.isDirectory())
		return false;

	if (!EnsureDirectoryStructureForDataImages(K_HEIGHTMAP_FOLDER_NAME))
		return false;
	if (!EnsureDirectoryStructureForDataImages(K_WATERMAP_FOLDER_NAME))
		return false;
	if (!EnsureDirectoryStructureForDataImages(K_RAINFLOWMAP_FOLDER_NAME))
		return false;
	if (!EnsureDirectoryStructureForDataImages(K_SEDIMENTMAP_FOLDER_NAME))
		return false;

	return true;
}
public String GetDirectory()
{
	String dir = parent.GetDirectory() + File.separator;
	dir += WorldMap.K_REGIONS_FOLDER_NAME + File.separator;
	dir += GetDirectoryName();
	return dir;
}
public String GetDirectoryName()
{
	return "RegionalMap_" + GetWorldX() + "_" + GetWorldY();
}*/
void 
RegionalMap::ScrollOriginOffsetForOptimalCoastliness()
{
	int bestK = 0;
	int numLandTilesInBest = 0;
	for (int k = 0; k < 200; k++)
	{
		int modifiedOffset = ORIGIN_OFFSET + k;
		int numTiles = 0;

		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				double x = 1.0 * i / DIMENSION + modifiedOffset;
				double y = 1.0 * j / DIMENSION + modifiedOffset;
				if (Perlin::oceans.UnderThreshold(x, y))
					numTiles++;
			}
		}
		if (numTiles > numLandTilesInBest)
		{
			bestK = k;
			numLandTilesInBest = numTiles;
		}
	}
	ORIGIN_OFFSET = ORIGIN_OFFSET + bestK;
}
SamplePoint* 
RegionalMap::GetSamplePoint(int index)
{
	if (index < 0 || index >= voronoiList.size())
		return nullptr;
	return voronoiList[index];
}
void 
RegionalMap::CalculateAllVoronoiAdjacencies()
{
	for (int i = 0; i < voronoiList.size(); i++)
	{
		voronoiList[i]->CalculateAdjacencies();
	}
}
bool 
RegionalMap::IsRidge(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	SamplePoint* vp = GetNearest(x, y);
	if (vp == nullptr)
		return false;

	array<SamplePoint*, 3> triangle = VoronoiAlgorithms::FindContainingSampleTriangle(x, y, vp);
	if (triangle[0] == nullptr)
		return false;

	for (SamplePoint* a : triangle)
	{
		for (SamplePoint* b : triangle)
		{
			if (a == b)
				continue;
			MeshConnection* m = a->GetConnection(b);
			if (m == nullptr)
				continue;
			if (!m->IsRidgeline())
				continue;

			double pToAX = (a->x - x);
			double pToAY = (a->y - y);
			double aToBX = (b->x - a->x);
			double aToBY = (b->y - a->y);
			double pToBX = (b->x - x);
			double pToBY = (b->y - y);
			double aDist = sqrt(pToAX * pToAX + pToAY * pToAY);
			double bDist = sqrt(pToBX * pToBX + pToBY * pToBY);

			double perpX = aToBY;
			double perpY = -aToBX;
			perpX /= m->GetLength();
			perpY /= m->GetLength();
			double perpDist = abs(perpX * pToAX + perpY * pToAY);

			if (aDist < 0.0001 || bDist < 0.0001 || perpDist < 0.0001)
				return true;
		}
	}
	return false;
}
double 
RegionalMap::DrainageBasedRiverPercent(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	SamplePoint* vp = GetNearest(x, y);
	if (vp == nullptr)
		return 0;

	array<MeshPoint*, 3> triangle = VoronoiAlgorithms::FindContainingTriangle(x, y, vp);
	if (triangle[0] == nullptr)
		return 0;
	double rain = 10000 * 4;
	double highestPercent = 0;
	for (MeshPoint* v : triangle)
	{
		double pToVX = (v->x - x);
		double pToVY = (v->y - y);
		double pToV = sqrt(pToVX * pToVX + pToVY * pToVY);
		highestPercent = max(highestPercent, RiverFlowPercent(rain * v->GetDrainageArea(), pToV));
		if (v->GetDrainageSink() != nullptr)
		{
			double vToOX = v->GetDrainageSink()->x - v->x;
			double vToOY = v->GetDrainageSink()->y - v->y;
			if (vToOX * pToVX + vToOY * pToVY < 0)
			{
				double perpX = vToOY;
				double perpY = -vToOX;
				perpX /= v->DistTo(v->GetDrainageSink());
				perpY /= v->DistTo(v->GetDrainageSink());
				double perpDist = abs(perpX * pToVX + perpY * pToVY);
				highestPercent = max(highestPercent, RiverFlowPercent(rain * v->GetDrainageArea(), perpDist));
			}
		}
		for (MeshPoint* s : v->GetDrainageSources())
		{
			double pToSX = (s->x - x);
			double pToSY = (s->y - y);
			double pToS = sqrt(pToSX * pToSX + pToSY * pToSY);
			highestPercent = max(highestPercent, RiverFlowPercent(rain * s->GetDrainageArea(), pToS));

			double sToOX = v->x - s->x;
			double sToOY = v->y - s->y;
			if (sToOX * pToVX + sToOY * pToVY > 0)
			{
				double perpX = sToOY;
				double perpY = -sToOX;
				perpX /= s->DistTo(v);
				perpY /= s->DistTo(v);
				double perpDist = abs(perpX * pToSX + perpY * pToSY);
				highestPercent = max(highestPercent, RiverFlowPercent(rain * s->GetDrainageArea(), perpDist));
			}
		}
	}
	return highestPercent;
}
double 
RegionalMap::RiverPercent(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	SamplePoint* vp = GetNearest(x, y);
	if (vp == nullptr)
		return 0;

	array<SamplePoint*, 3> triangle = VoronoiAlgorithms::FindContainingSampleTriangle(x, y, vp);
	if (triangle[0] == nullptr)
		return 0;

	double highestPercent = 0;
	for (SamplePoint* v : triangle)
	{
		double pToVX = (v->x - x);
		double pToVY = (v->y - y);
		double pToV = sqrt(pToVX * pToVX + pToVY * pToVY);
		highestPercent = max(highestPercent, RiverFlowPercent(v->GetRiverFlow(), pToV));
		if (v->GetRiverOutlet() != nullptr)
		{
			double vToOX = v->GetRiverOutlet()->x - v->x;
			double vToOY = v->GetRiverOutlet()->y - v->y;
			if (vToOX * pToVX + vToOY * pToVY < 0)
			{
				double perpX = vToOY;
				double perpY = -vToOX;
				perpX /= v->DistTo(v->GetRiverOutlet());
				perpY /= v->DistTo(v->GetRiverOutlet());
				double perpDist = abs(perpX * pToVX + perpY * pToVY);
				highestPercent = max(highestPercent, RiverFlowPercent(v->GetRiverFlow(), perpDist));
			}
		}
		for (SamplePoint* s : v->GetRiverInlets())
		{
			double pToSX = (s->x - x);
			double pToSY = (s->y - y);
			double pToS = sqrt(pToSX * pToSX + pToSY * pToSY);
			highestPercent = max(highestPercent, RiverFlowPercent(s->GetRiverFlow(), pToS));

			double sToOX = v->x - s->x;
			double sToOY = v->y - s->y;
			if (sToOX * pToVX + sToOY * pToVY > 0)
			{
				double perpX = sToOY;
				double perpY = -sToOX;
				perpX /= s->DistTo(v);
				perpY /= s->DistTo(v);
				double perpDist = abs(perpX * pToSX + perpY * pToSY);
				highestPercent = max(highestPercent, RiverFlowPercent(s->GetRiverFlow(), perpDist));
			}
		}
	}
	return highestPercent;
}
double 
RegionalMap::RiverFlowPercent(double flow, double distanceFromFlow)
{
	double invisibleFlowThreshold = 20;//20;
	double fullVisibleFlowThreshold = 800;//800;
	if (flow < invisibleFlowThreshold)
		return 0;
	if (distanceFromFlow > MIN_VORONOI_DIST)
		return 0;
	double distScaled = distanceFromFlow / 0.0002;
	double flowWidth = pow(flow, 0.5);
	double zeroPercentDist = flowWidth - pow(invisibleFlowThreshold, 0.5);
	double fullPercentDist = flowWidth - pow(fullVisibleFlowThreshold, 0.5);
	if (distScaled > zeroPercentDist)
		return 0;
	if (distScaled < fullPercentDist)
		return 1;
	return (zeroPercentDist - distScaled) / (zeroPercentDist - fullPercentDist);
}

unordered_map<LocalMap*, bool> 
RegionalMap::PrepareForExtensiveEditing()
{
	unordered_map<LocalMap*, bool> used;
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			cout << "Preparing Local Map " << i << ", " << j << endl;
			LocalMap* lm = topography[i * DIMENSION + j];
			bool alreadyActive = lm->PrepareForEditing(true, true, true);
			used[lm] = alreadyActive;
		}
	}
	cout << endl;
	for (int i = 0; i < DIMENSION; i++)
	{
		LocalMap* nEdge = topography[i * DIMENSION];
		LocalMap* sEdge = topography[i * DIMENSION + DIMENSION - 1];
		LocalMap* north = nEdge->GetNorth();
		LocalMap* south = sEdge->GetSouth();
		if (north != nullptr)
		{
			bool alreadyActive = north->PrepareForEditing(true, true, true);
			used[north] = alreadyActive;
		}
		if (south != nullptr)
		{
			bool alreadyActive = south->PrepareForEditing(true, true, true);
			used[south] = alreadyActive;
		}
		if (i == 0 && north != nullptr)
		{
			LocalMap* northWest = north->GetWest();
			if (northWest != nullptr)
			{
				bool alreadyActive = northWest->PrepareForEditing(true, true, true);
				used[northWest] = alreadyActive;
			}
		}
		if (i == 0 && south != nullptr)
		{
			LocalMap* southWest = south->GetWest();
			if (southWest != nullptr)
			{
				bool alreadyActive = southWest->PrepareForEditing(true, true, true);
				used[southWest] = alreadyActive;
			}
		}
		if (i == DIMENSION - 1 && north != nullptr)
		{
			LocalMap* northEast = north->GetEast();
			if (northEast != nullptr)
			{
				bool alreadyActive = northEast->PrepareForEditing(true, true, true);
				used[northEast] = alreadyActive;
			}
		}
		if (i == DIMENSION - 1 && south != nullptr)
		{
			LocalMap* southEast = south->GetEast();
			if (southEast != nullptr)
			{
				bool alreadyActive = southEast->PrepareForEditing(true, true, true);
				used[southEast] = alreadyActive;
			}
		}
	}
	for (int j = 0; j < DIMENSION; j++)
	{
		LocalMap* wEdge = topography[j];
		LocalMap* eEdge = topography[(DIMENSION - 1) * DIMENSION + j];
		LocalMap* west = wEdge->GetWest();
		LocalMap* east = eEdge->GetEast();
		if (west != nullptr)
		{
			bool alreadyActive = west->PrepareForEditing(true, true, true);
			used[west] = alreadyActive;
		}
		if (east != nullptr)
		{
			bool alreadyActive = east->PrepareForEditing(true, true, true);
			used[east] = alreadyActive;
		}
	}
	return used;
}
void 
RegionalMap::RunFullPhasedErosion()
{

	vector<LocalMap*> targets;
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			LocalMap* lm = topography[i * DIMENSION + j];
			targets.push_back(lm);
		}
	}
	unordered_map<LocalMap*, bool> used = PrepareForExtensiveEditing();
	cout << "Preparing Pixels" << endl;
	unordered_map<LocalMap*, bool> usedForPixels;
	//vector<LocalMap::Pixel> allPixels =
	//	LocalTerrainAlgorithms::BeginPixelOperation(usedForPixels, targets, true, true, true, true, true);

	for (int n = 0; n < 20; n++)
	{
		cout << "*********" << endl;
		cout << "LAPLACE IT " << n << endl;
		cout << "*********" << endl << endl;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				int blockSize = DIMENSION / 2;
				for(int index = 0; index < blockSize * blockSize; index++)
				{
					int myI = index / blockSize;
					int myJ = index % blockSize;
					int realI = i + myI * 2;
					int realJ = j + myJ * 2;
					LocalMap* target = topography[realI * DIMENSION + realJ];
					target->LaplacianErosionIteration(1);
				}
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << "RAIN 1 BLOCK " << i << ", " << j << endl;
			cout << "********" << endl << endl;
			int blockSize = DIMENSION / 4;
			for(int index = 0; index < blockSize * blockSize; index++)
			{
				int myI = index / blockSize;
				int myJ = index % blockSize;
				int realI = i + myI * 4;
				int realJ = j + myJ * 4;
				LocalMap* target = topography[realI * DIMENSION + realJ];
				for (int drop = 0; drop < 10000; drop++)
				{
					WaterDroplet nova(target, 1, used, false);
					bool okay = true;
					while (okay)
					{
						okay = nova.OneErosionStep();
					}
				}
			}
		}
	}
	cout << "********" << endl;
	cout << "FLUVIAL THERMAL 1" << endl;
	cout << "********" << endl << endl;
	//LocalTerrainAlgorithms::ThermalFluvialErosion(allPixels, 10);

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << "********" << endl;
			cout << "RAIN 2 BLOCK " << i << ", " << j << endl;
			cout << "********" << endl << endl;
			int blockSize = DIMENSION / 4;

			for(int index = 0; index < blockSize * blockSize; index++)
			{
				int myI = index / blockSize;
				int myJ = index % blockSize;
				int realI = i + myI * 4;
				int realJ = j + myJ * 4;
				LocalMap* target = topography[realI * DIMENSION + realJ];
				for (int drop = 0; drop < 10000; drop++)
				{
					WaterDroplet nova(target, 1, used, false);
					bool okay = true;
					while (okay)
					{
						okay = nova.OneErosionStep();
					}
				}
			}
		}
	}

	cout << "********" << endl;
	cout << "FINAL HYDROLOGY" << endl;
	cout << "********" << endl << endl;
	//LocalTerrainAlgorithms::GuaranteeConsistentHydrology(allPixels);

	//LocalTerrainAlgorithms::EndPixelOperation(usedForPixels, true, true, true, true, true);

	for (auto lm : used)
	{
		lm.first->CompleteEditing(true, true, true, !lm.second);
	}
}

void 
RegionalMap::RunFullLaplacianErosion()
{
	int numIterations = 50;

	unordered_map<LocalMap*, bool> used = PrepareForExtensiveEditing();

	for (int n = 0; n < numIterations; n++)
	{
		cout << "********" << endl;
		cout << "LPLIT " << n << endl;
		cout << "********" << endl;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				cout << endl << "Starting Laplace Erosion on block " << i << ", " << j << endl << endl;
				int blockSize = DIMENSION / 2;

				for(int index = 0; index < blockSize * blockSize; index++)
				{
					int myI = index / blockSize;
					int myJ = index % blockSize;
					int realI = i + myI * 2;
					int realJ = j + myJ * 2;
					LocalMap* target = topography[realI * DIMENSION + realJ];
					target->LaplacianErosionIteration(1);
				}
			}
		}
	}

	for (auto lm : used)
	{
		lm.first->CompleteEditing(true, true, true, !lm.second);
	}
}
void 
RegionalMap::RunFullRain()
{
	int numDroplets = 5000;

	unordered_map<LocalMap*, bool> used = PrepareForExtensiveEditing();

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << endl << endl << "Starting rain on block " << i << ", " << j << endl << endl << endl;
			int blockSize = DIMENSION / 4;

			for (int index = 0; index < blockSize * blockSize; index++)
			{
				int myI = index / blockSize;
				int myJ = index % blockSize;
				int realI = i + myI * 4;
				int realJ = j + myJ * 4;
				cout << "Raining on map " << realI << ", " << realJ << endl;
				LocalMap* target = topography[realI * DIMENSION + realJ];
				for (int drop = 0; drop < numDroplets; drop++)
				{
					WaterDroplet nova(target, 1, used, false);
					bool okay = true;
					while (okay)
					{
						okay = nova.OneErosionStep();
					}
				}
			}
		}
	}
	for (auto lm : used)
	{
		lm.first->CompleteEditing(true, true, true, !lm.second);
	}
}
double 
RegionalMap::GetElevation(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	LocalMap::Coordinate lmc = GetLocalMapAt(x - GetWorldX(), y - GetWorldY());
	if (lmc.src == nullptr)
		return 0;
	return lmc.src->GetHeight(lmc.x, lmc.y);
}

WatermapValue 
RegionalMap::IsWaterPoint(double wX, double wY)
{
	double xub = wX / RegionalMap::DIMENSION;
	double yub = wY / RegionalMap::DIMENSION;
	double xb = xub;
	double yb = yub;
	if (Switches::USE_BLURRING)
	{
		double blurredX = Perlin::blurX.Get(xub, yub);
		double blurredY = Perlin::blurY.Get(xub, yub);
		xb += Perlin::BLUR_DISTANCE * blurredX / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
		yb += Perlin::BLUR_DISTANCE * blurredY / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
	}

	SamplePoint* vp = GetNearest(xb, yb);
	if (vp == nullptr)
		return WatermapValue::Unknown;

	array<MeshPoint*, 3> triangle = VoronoiAlgorithms::FindContainingTriangle(xb, yb, vp);
	if (triangle[0] == nullptr)
	{
		if (vp->IsOcean())
			return WatermapValue::Ocean;
		else if (vp->IsInlandLake())
			return WatermapValue::Lake;
		else
			return WatermapValue::Unknown;
	}
	bool anyOcean = false;
	bool anyLake = false;
	bool anyLand = false;
	for (MeshPoint* mp : triangle)
	{
		if (mp->IsOcean())
			anyOcean = true;
		else if (mp->IsInlandLake())
			anyLake = true;
		else
			anyLand = true;
	}
	if (anyLand)
		return WatermapValue::NotWater;
	else if (anyLake)
		return WatermapValue::Lake;
	else if (anyOcean)
		return WatermapValue::Ocean;
	else
		return WatermapValue::Unknown;
}

//this is a complicated bit of math which comes from a mathematical machine
//I've sketched the operation of that machine on desmos; see this link:
//https://www.desmos.com/calculator/f6e0udqb1c
//The idea is that the sediment is going to cut steps into the terrain
//so that when it erodes away, a tiered, staircase structure remains.
//So step height controls the height of those steps;
//but these steps aren't perfectly vertical: they are divided into a flat
//section and a steep section.
//flatLength is a number between 0-1, and controls how long the flat section is.
//A good number would be 0.8-0.9.
//flatness controls how flat that section is; a value of 0 means that this
//staircase structure won't exist at all, a value of 1 makes the flat section
//completely, 100% flat. I find a value of 0.9-0.95 reasonable.
//Finally, the steep section will itself be divided into substeps
//numSubsteps controls how many. This number must be at least 1; a higher number
//actually makes things smoother. 2-4 gives a series of blocky substeps.
double 
RegionalMap::SedimentStairCalculation(
	double terrainHeight, 
	double stepHeight, 
	double flatLength, 
	double flatness, 
	int numSubsteps)
{
	double sedimentInStep = fmod(terrainHeight, stepHeight);
	while (sedimentInStep < 0)
		sedimentInStep += stepHeight;

	if (sedimentInStep <= flatLength * stepHeight)
		return sedimentInStep * flatness;

	double substepHeight = stepHeight - flatLength * stepHeight;
	substepHeight /= numSubsteps;

	double substep = fmod((sedimentInStep - flatLength * stepHeight), substepHeight);
	while (substep < 0)
		substep += substepHeight;

	double substepModifier = flatLength * flatness / (1 - flatLength);
	substepModifier += 1;
	substepModifier *= substep;

	double steepHeight = flatness * flatLength * (stepHeight - sedimentInStep);
	steepHeight /= (1 - flatLength);

	return steepHeight + substepModifier;
}
double 
RegionalMap::CalculateBaseSedimentDepth(double wX, double wY)
{
	double xub = wX / RegionalMap::DIMENSION;
	double yub = wY / RegionalMap::DIMENSION;
	double xb = xub;
	double yb = yub;
	if (Switches::USE_BLURRING)
	{
		double blurredX = Perlin::blurX.Get(xub, yub);
		double blurredY = Perlin::blurY.Get(xub, yub);
		xb += Perlin::BLUR_DISTANCE * blurredX / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
		yb += Perlin::BLUR_DISTANCE * blurredY / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
	}

	SamplePoint* vp = GetNearest(xb, yb);
	if (vp == nullptr)
		return 0;

	double elev = 0;
	array<MeshPoint*, 3> meshTriangle = VoronoiAlgorithms::FindContainingTriangle(xb, yb, vp);
	if (meshTriangle[0] == nullptr)
		elev = vp->GetElevation();
	else
	{
		glm::dvec3 lerp = VoronoiAlgorithms::BarycentricCoordinates(xb, yb, meshTriangle);
		for (int i = 0; i < 3; i++)
		{
			double elevCont = meshTriangle[i]->GetElevation();
			if (elevCont == numeric_limits<int>::max())
				elevCont = 0;
			elev += elevCont * lerp[i];
		}
	}

	double stepWidth = Perlin::sedimentStepDelta.Get(xub, yub);
	double signum = (stepWidth > 0.0) ? 1.0 : -1.0;
	stepWidth = abs(stepWidth);
	stepWidth = pow(stepWidth, 0.5);
	stepWidth *= signum;
	stepWidth += 0.8;
	if (stepWidth < 0)
		stepWidth = 0;
	stepWidth *= 700;
	stepWidth += 40;

	double sedDep = SedimentStairCalculation(elev, stepWidth, 0.65, 0.8, 3);

	double mask = Perlin::sedimentStepMask.Get(xub, yub);
	if (mask > -0.17)
		mask = (mask + 0.17) * 500;
	else
		mask = 0;

	//sedDep = Math.min(sedDep, mask);

	array<SamplePoint*, 3> triangle = VoronoiAlgorithms::FindContainingSampleTriangle(xb, yb, vp);
	if (triangle[0] == nullptr)
	{
		return vp->GetBaseSedimentDepth() + sedDep;
	}

	glm::dvec3 lerp = VoronoiAlgorithms::BarycentricCoordinates(xb, yb, triangle);
	double lerpedVal = 0;
	for (int i = 0; i < 3; i++)
	{
		lerpedVal += triangle[i]->GetBaseSedimentDepth() * lerp[i];
	}
	return lerpedVal + sedDep;
}
double 
RegionalMap::CalculateElevation(double wX, double wY)
{
	double xub = wX / RegionalMap::DIMENSION;
	double yub = wY / RegionalMap::DIMENSION;
	double xb = xub;
	double yb = yub;
	if (Switches::USE_BLURRING)
	{
		double blurredX = Perlin::blurX.Get(xub, yub);
		double blurredY = Perlin::blurY.Get(xub, yub);
		xb += Perlin::BLUR_DISTANCE * blurredX / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
		yb += Perlin::BLUR_DISTANCE * blurredY / (RegionalMap::DIMENSION * LocalMap::METER_DIM);
	}

	SamplePoint* vp = GetNearest(xb, yb);
	if (vp == nullptr)
		return 0;

	array<MeshPoint*, 3> triangle = VoronoiAlgorithms::FindContainingTriangle(xb, yb, vp);
	if (triangle[0] == nullptr)
		return vp->GetElevation();

	double elev = 0;
	glm::dvec3 lerp = VoronoiAlgorithms::BarycentricCoordinates(xb, yb, triangle);
	double perlinContribs[Perlin::elevDeltaLength];
	for (int i = 0; i < 3; i++)
	{
		double elevCont = triangle[i]->GetElevation();
		vector<double> deltaContribs = triangle[i]->GetPerlinElevDiffs();
		if (elevCont == numeric_limits<int>::max())
			elevCont = 0;
		elev += elevCont * lerp[i];
		for (int j = 0; j < Perlin::elevDeltaLength; j++)
			perlinContribs[j] += lerp[i] * deltaContribs[j];
	}

	for (int i = 0; i < Perlin::elevDeltaLength; i++)
	{
		double del = Perlin::elevDeltas[i].Get(xub, yub);
		double ctr = perlinContribs[i];
		double amp = Perlin::elevDeltaScales[i];
		elev += amp * ctr * del;
	}
	if (elev < 0)
		return 0;

	WatermapValue isWater = IsWaterPoint(wX, wY);
	if (isWater != WatermapValue::NotWater)
		return elev;
	double rockJitter = Perlin::rockyJitters.Get(xub, yub);
	rockJitter *= perlinContribs[1];
	rockJitter *= rockJitter;
	if (elev + rockJitter * Perlin::rockJitterScale < 0)
		return 0;
	return elev + rockJitter * Perlin::rockJitterScale;
}
double 
RegionalMap::GetElevationLaplacian(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	LocalMap::Coordinate lmc = GetLocalMapAt(x - GetWorldX(), y - GetWorldY());
	if (lmc.src == nullptr)
		return 0;
	return lmc.src->GetHeightLaplacian(lmc.x, lmc.y);
}
glm::dvec2 
RegionalMap::GetElevationGradient(double wX, double wY)
{
	double x = wX / RegionalMap::DIMENSION;
	double y = wY / RegionalMap::DIMENSION;
	LocalMap::Coordinate lmc = GetLocalMapAt(x - GetWorldX(), y - GetWorldY());
	if (lmc.src == nullptr)
		return glm::dvec2(0, 0);
	return lmc.src->GetHeightGradient(lmc.x, lmc.y);
}

vector<SamplePoint*> 
RegionalMap::GetNearestN(double wX, double wY, int n)
{
	vector<SamplePoint*> confirmed;
	int i = (int)((wX - GetWorldX()) * VORONOI_DIM);
	int j = (int)((wY - GetWorldY()) * VORONOI_DIM);
	map<double, SamplePoint*> found;
	int d = 1;
	while (true)
	{
		for (int dI = i - d; dI <= i + d; dI++)
		{
			SamplePoint* comp = GetCell(dI, j - d);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			found[dX * dX + dY * dY] = comp;
		}
		for (int dI = i - d; dI <= i + d; dI++)
		{
			SamplePoint* comp = GetCell(dI, j + d);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			found[dX * dX + dY * dY] = comp;
		}
		for (int dJ = j - d; dJ <= j + d; dJ++)
		{
			SamplePoint* comp = GetCell(i - d, dJ);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			found[dX * dX + dY * dY] = comp;
		}
		for (int dJ = j - d; dJ <= j + d; dJ++)
		{
			SamplePoint* comp = GetCell(i + d, dJ);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			found[dX * dX + dY * dY] = comp;
		}
		//min next ring distance possible
		double mnrd = (d * 1.0) / VORONOI_DIM;
		while (!found.empty() && found.begin()->first < mnrd * mnrd)
		{
			SamplePoint* value = found.begin()->second;
			found.erase(found.begin());
			confirmed.push_back(value);
			if (confirmed.size() >= n)
				break;
		}
		if (confirmed.size() >= n)
			break;
		d++;
	}
	return confirmed;
}

SamplePoint* 
RegionalMap::GetNearest(double wX, double wY)
{
	int i = (int)((wX - GetWorldX()) * VORONOI_DIM);
	int j = (int)((wY - GetWorldY()) * VORONOI_DIM);

	SamplePoint* best = GetCell(i, j);
	double bestDistSqr = 0;
	if (best != nullptr)
	{
		double dX = wX - best->x;
		double dY = wY - best->y;
		bestDistSqr = dX * dX + dY * dY;
	}
	int d = 1;
	while (true)
	{
		for (int dI = i - d; dI <= i + d; dI++)
		{
			SamplePoint* comp = GetCell(dI, j - d);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			if (best == nullptr || bestDistSqr > dX * dX + dY * dY)
			{
				best = comp;
				bestDistSqr = dX * dX + dY * dY;
			}
		}
		for (int dI = i - d; dI <= i + d; dI++)
		{
			SamplePoint* comp = GetCell(dI, j + d);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			if (best == nullptr || bestDistSqr > dX * dX + dY * dY)
			{
				best = comp;
				bestDistSqr = dX * dX + dY * dY;
			}
		}
		for (int dJ = j - d; dJ <= j + d; dJ++)
		{
			SamplePoint* comp = GetCell(i - d, dJ);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			if (best == nullptr || bestDistSqr > dX * dX + dY * dY)
			{
				best = comp;
				bestDistSqr = dX * dX + dY * dY;
			}
		}
		for (int dJ = j - d; dJ <= j + d; dJ++)
		{
			SamplePoint* comp = GetCell(i + d, dJ);
			if (comp == nullptr)
				continue;
			double dX = wX - comp->x;
			double dY = wY - comp->y;
			if (best == nullptr || bestDistSqr > dX * dX + dY * dY)
			{
				best = comp;
				bestDistSqr = dX * dX + dY * dY;
			}
		}
		//min next ring distance possible
		double mnrd = (d * 1.0) / VORONOI_DIM;
		if (best != nullptr && bestDistSqr < mnrd * mnrd)
			break;
		if (d > 5)
			break;
		d++;
	}

	return best;
}
bool 
RegionalMap::ExistingPointNearby(SamplePoint* cand)
{
	int i = (int)((cand->x - GetWorldX()) * VORONOI_DIM);
	int j = (int)((cand->y - GetWorldY()) * VORONOI_DIM);
	for (int dI = -2; dI <= 2; dI++)
		for (int dJ = -2; dJ <= 2; dJ++)
		{
			SamplePoint* comp = GetCell(i + dI, j + dJ);
			if (comp == nullptr)
				continue;
			double dX = cand->x - comp->x;
			double dY = cand->y - comp->y;
			if (dX * dX + dY * dY < MIN_VORONOI_DIST * MIN_VORONOI_DIST)
				return true;
		}
	return false;
}
SamplePoint* 
RegionalMap::GetCell(int i, int j)
{
	if (i < 0)
	{
		RegionalMap* alt = GetWest();
		if (alt == nullptr)
			return nullptr;
		return alt->GetCell(i + VORONOI_DIM, j);
	}
	if (j < 0)
	{
		RegionalMap* alt = GetNorth();
		if (alt == nullptr)
			return nullptr;
		return alt->GetCell(i, j + VORONOI_DIM);
	}
	if (i >= VORONOI_DIM)
	{
		RegionalMap* alt = GetEast();
		if (alt == nullptr)
			return nullptr;
		return alt->GetCell(i - VORONOI_DIM, j);
	}
	if (j >= VORONOI_DIM)
	{
		RegionalMap* alt = GetSouth();
		if (alt == nullptr)
			return nullptr;
		return alt->GetCell(i, j - VORONOI_DIM);
	}
	return terrainCells[i * VORONOI_DIM + j];
}
bool 
RegionalMap::SetPoint(SamplePoint* p)
{
	if (p->x - GetWorldX() < 0 || p->y - GetWorldY() < 0)
		return false;
	int i = (int)((p->x - GetWorldX()) * VORONOI_DIM);
	int j = (int)((p->y - GetWorldY()) * VORONOI_DIM);
	if (i >= VORONOI_DIM || j >= VORONOI_DIM)
		return false;
	if (terrainCells[i * VORONOI_DIM + j] == nullptr)
	{
		terrainCells[i * VORONOI_DIM + j] = p;
		int index = voronoiList.size();
		p->SetContainerIndex(index);
		voronoiList.push_back(p);
		p->AssignVoronoiTerrainType();
		return true;
	}
	return false;
}

/*public void Render(double d, Graphics2D g2)
{
	if (!readyToRender)
	{
		g2.setColor(new Color(50, 50, 50));
		g2.fillRect(0, 0, (int)d, (int)d);
		return;
	}
	AffineTransform saved = g2.getTransform();
	double tileD = d / DIMENSION;
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			if (g2.hitClip(0, 0, (int)tileD, (int)tileD))
			{
				LocalMap lm = topography[i * DIMENSION + j];
				lm.Render((int)tileD, g2);
				if (Switches.OUTLINE_MAPS)
				{
					g2.setColor(Color.red);
					g2.drawRect(0, 0, (int)tileD, (int)tileD);
				}
			}
			g2.translate(0, tileD);

		}
		g2.translate(tileD, -1 * DIMENSION * tileD);
	}
	g2.setTransform(saved);
	if (Switches.OUTLINE_MAPS)
	{
		g2.setColor(Color.green);
		g2.drawRect(0, 0, (int)d, (int)d);
	}


	if (Switches.PAINT_VORONOI_CENTERS)
	{
		double mvd = MIN_VORONOI_DIST * d;
		int dotMvd = 6;
		if (mvd > 12)
		{
			for (int i = 0; i < voronoiList.size(); i++)
			{
				SamplePoint p = voronoiList.get(i);
				double x = p.x - GetWorldX();
				double y = p.y - GetWorldY();
				x *= d;
				y *= d;
				int xL = (int)(x - mvd * 0.5);
				int yL = (int)(y - mvd * 0.5);
				if (!g2.hitClip(xL, yL, (int)mvd, (int)mvd))
					continue;

				if (mvd > 24)
				{
					g2.setColor(Color.orange);
					g2.drawOval(xL, yL, (int)(mvd), (int)(mvd));
				}

				xL = (int)(x - dotMvd * 0.5);
				yL = (int)(y - dotMvd * 0.5);
				if (p.IsWaterPoint())
				{
					g2.setColor(Color.blue);
					g2.fillOval(xL - 1, yL - 1, (int)(dotMvd)+2, (int)(dotMvd)+2);
				}
				g2.setColor(Color.black);
				g2.fillOval(xL, yL, (int)(dotMvd), (int)(dotMvd));
				if (mvd < 48)
					continue;
				for (SamplePoint s : p.GetAdjacentSamples())
				{
					MeshConnection con = p.GetConnection(s);
					MeshPoint m = con.GetMid();
					if (m == null)
						continue;
					x = m.x - GetWorldX();
					y = m.y - GetWorldY();
					x *= d;
					y *= d;

					xL = (int)(x - dotMvd * 0.5);
					yL = (int)(y - dotMvd * 0.5);
					if (!g2.hitClip(xL, yL, (int)(dotMvd), (int)(dotMvd)))
						continue;

					if (m.IsWaterPoint())
					{
						g2.setColor(Color.blue);
						g2.fillOval(xL - 1, yL - 1, (int)(dotMvd)+2, (int)(dotMvd)+2);
					}
					g2.setColor(Color.green);
					g2.fillOval(xL, yL, (int)(dotMvd), (int)(dotMvd));

					if (mvd < 96)
						continue;
					for (MeshPoint mAdj : m.GetAdjacent())
					{
						MeshConnection mCon = m.GetConnection(mAdj);
						if (mCon == null)
							continue;
						MeshPoint mm = mCon.GetMid();
						if (mm == null)
							continue;
						x = mm.x - GetWorldX();
						y = mm.y - GetWorldY();
						x *= d;
						y *= d;

						xL = (int)(x - dotMvd * 0.5);
						yL = (int)(y - dotMvd * 0.5);
						if (!g2.hitClip(xL, yL, (int)(dotMvd), (int)(dotMvd)))
							continue;

						if (mm.IsWaterPoint())
						{
							g2.setColor(Color.blue);
							g2.fillOval(xL - 1, yL - 1, (int)(dotMvd)+2, (int)(dotMvd)+2);
						}
						g2.setColor(Color.red);
						g2.fillOval(xL, yL, (int)(dotMvd), (int)(dotMvd));

						if (mvd < 128)
							continue;
						for (MeshPoint mmadj : mm.GetAdjacent())
						{
							MeshConnection mmCon = mm.GetConnection(mmadj);
							if (mmCon == null)
								continue;
							MeshPoint mmm = mmCon.GetMid();
							if (mmm == null)
								continue;
							x = mmm.x - GetWorldX();
							y = mmm.y - GetWorldY();
							x *= d;
							y *= d;

							xL = (int)(x - dotMvd * 0.5);
							yL = (int)(y - dotMvd * 0.5);
							if (!g2.hitClip(xL, yL, (int)(dotMvd), (int)(dotMvd)))
								continue;
							g2.setColor(Color.blue);
							g2.fillOval(xL, yL, (int)(dotMvd), (int)(dotMvd));
						}
					}
				}

			}
		}

	}
}*/
LocalMap* 
RegionalMap::GetLocalMapAt(int x, int y)
{
	if (x < 0 || y < 0)
		return nullptr;
	if (x > DIMENSION || y > DIMENSION)
		return nullptr;
	return topography[x * DIMENSION + y];
}
RegionalMap*
RegionalMap::GetNorth()
{
	return parent->GetRegion(x, y - 1);
}
RegionalMap*
RegionalMap::GetEast()
{
	return parent->GetRegion(x + 1, y);
}
RegionalMap*
RegionalMap::GetWest()
{
	return parent->GetRegion(x - 1, y);
}
RegionalMap*
RegionalMap::GetSouth()
{
	return parent->GetRegion(x, y + 1);
}

//The assumption here is that the x, y are in the [0, 1] interval
LocalMap::Coordinate 
RegionalMap::GetLocalMapAt(double x, double y)
{
	x *= DIMENSION;
	y *= DIMENSION;
	int xIndex = (int)x;
	int yIndex = (int)y;
	if (xIndex < 0 || xIndex >= DIMENSION)
		return LocalMap::Coordinate(0, 0, nullptr);
	if (yIndex < 0 || yIndex >= DIMENSION)
		return LocalMap::Coordinate(0, 0, nullptr);

	x -= xIndex;
	y -= yIndex;

	return LocalMap::Coordinate(x, y, topography[xIndex * DIMENSION + yIndex]);
}
