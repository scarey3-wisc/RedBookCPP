#pragma once
#include <shared_mutex>
#include <unordered_map>
#include <list>
#include <array>
#include "DataImageSize.h"
#include <string>
#include <iostream>

template<typename T>
class LinkedHashSet
{
public:

	//------------------------------------------------------------------------------

	void Clear()
	{
		map.clear();
		list.clear();
	}

	//------------------------------------------------------------------------------

	T& KickLast()
	{
		auto last = list.back();
		list.pop_back();
		map.erase(last);
		return last;
	}

	//------------------------------------------------------------------------------

	int GetSize() const
	{
		return map.size();
	}

	//------------------------------------------------------------------------------

	void Remove(const T& value)
	{
		auto it = map.find(value);
		if (it == map.end())
			return;
		list.erase(it->second);
		map.erase(it);
	}

	//------------------------------------------------------------------------------

	void Add(const T& value)
	{
		auto it = map.find(value);
		if (it != map.end())
		{
			list.splice(list.begin(), list, it->second);
			return;
		}
		list.push_front(value);
		map[value] = list.begin();
	}

	//------------------------------------------------------------------------------

	bool Contains(const T& value, bool refresh = true) const
	{
		if (refresh)
		{
			auto it = map.find(value);
			if (it == map.end())
				return false;
			list.splice(list.begin(), list, it->second);
			return true;
		}
		return map.find(value) != map.end();
	}
private:
	std::unordered_map<T, typename std::list<T>::iterator> map;
	std::list<T> list;
};


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


template<typename Data, typename ImageType>
class DataImage;

template<typename A>
class DataProviderBase
{
public:
	virtual A GetNorth(); //meaning y < 0; we want to query the north image with y+1
	virtual A GetSouth(); //meaning y > 1; we want to query the south image with y-1
	virtual A GetWest(); //meaning x < 0
	virtual A GetEast(); //meaning x > 1
};


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


template<typename Data, typename ImageType>
class DataImageManager {
	friend class DataImage<Data, ImageType>;
public:
	DataImageManager(int capacities[NumDataImageDimensions])
	{
		for (int i = 0; i < NumDataImageDimensions; ++i)
		{
			this->capacities[i] = capacities[i];
		}
	}
	void SaveAll()
	{
		for (int i = 0; i < NumDataImageDimensions; i++)
		{
			//for (DataImage<B, A> di : caches.get(i))
			 //   di.SaveAllResolutions(false);
		}
	}
	int GetCurrentSize(int index)
	{
		if (index >= 0 && index < NumDataImageDimensions)
			return dataImages[index].GetSize();
		return 0;
	}
	int GetCapacity(int index)
	{
		if (index >= 0 && index < NumDataImageDimensions)
			return capacities[index];
		return 0;
	}

private:
	void FullKick(DataImage<Data, ImageType>* requester);

	void RequestCachePresenceAtCurrentResolution(DataImage<Data, ImageType>* requester);

	void KickingInsert(DataImage<Data, ImageType>* add, int index, std::vector<DataImage<Data, ImageType>*> kickRecord);


private:
	std::shared_mutex mutex;
	int capacities[NumDataImageDimensions] = {};
	std::array<LinkedHashSet<DataImage<Data, ImageType>*>, NumDataImageDimensions> dataImages;
};

template <typename Data, typename ImageType>
class DataImage
{
public:
	DataImage(
		const char* parentDir,
		const char* dirPrefix,
		const char* folderName,
		const char* fullName,
		DataProviderBase<ImageType>* source,
		DataImageManager<Data, ImageType>* mgr) :
		parentDir(parentDir), dirPrefix(dirPrefix), folderName(folderName), fullName(fullName),
		source(source), manager(mgr), changesRelativeToFile(false) {
	}



	std::string GetFullName() { return fullName; }

	double GetPixelSize()
	{
		if (GetCurrentResolution() == -1)
			return 0;
		int bestDim = DataImageDimensions[GetCurrentResolution()];
		if (bestDim == -1)
			return -1;
		double width = 1.0 / bestDim;
		return width;
	}
	void AnnounceChangesRelativeToFile() { changesRelativeToFile = true; }
	void RemoveFromManagement()
	{
		if (manager == nullptr)
			return;
		manager->FullKick(this);
		manager = nullptr;
	}
	void GiveToManagement(DataImageManager<Data, ImageType>* newManager)
	{
		//we're already under this management!
		if (manager == newManager)
			return;

		//We need to get out of the previous management
		if (manager != nullptr)
			RemoveFromManagement();

		newManager->RequestCachePresenceAtCurrentResolution(this);
		manager = newManager;
	}
	void ForceEditReady(bool verbose)
	{
		if (verbose)
			std::cout << "DataImage " << fullName << " getting ready to edit" << std::endl;
		DemandResolutionLevel(NumDataImageDimensions - 1);
	}

	virtual bool SaveImage(std::vector<Data>& image, std::string filename);



	bool SaveAllResolutions(bool force)
	{
		if (!force && !changesRelativeToFile)
			return true;
		int index = 0;
		for (std::vector<Data>& img : data)
		{
			std::string fileName = GetFileName(index);
			if (!SaveImage(img, fileName))
				return false;
			index++;
		}
		changesRelativeToFile = false;
		return true;
	}
	int GetCurrentResolution(){ return data.size() - 1; }
	void DemandResolution(int dim)
	{
		int targetResolution = 0;
		while (DataImageDimensions[targetResolution] < dim)
		{
			targetResolution++;
			if (targetResolution == NumDataImageDimensions)
			{
				targetResolution--;
				break;
			}
		}
		DemandResolutionLevel(targetResolution);
	}
	ImageType* GetNorth()
	{
		return source->GetNorth();
	}
	ImageType* GetSouth()
	{
		return source->GetSouth();
	}
	ImageType* GetWest()
	{
		return source->GetWest();
	}
	ImageType* GetEast()
	{
		return source->GetEast();
	}

protected:

	virtual std::vector<Data> Read(FILE* file);
	virtual std::vector<Data> Create(int dimIndex);

	void RemoveBestResolution()
	{
		data.erase(std::prev(data.end()));
	}

private:


	//res is an index in dimRange; so resolution 0 is likely 16x16 pixels
	void DemandResolutionLevel(int res)
	{
		//We're at or above that resolution, we're fine
		if (GetCurrentResolution() >= res)
			return;
		std::unique_lock lock(imgLock);
		while (GetCurrentResolution() < res)
		{
			int indexToAdd = GetCurrentResolution() + 1;
			//You've asked for a resolution level that we don't do; we're done here
			if (indexToAdd >= NumDataImageDimensions)
				break;

			//First, check if we have this next resolution on disk
			std::string fileName = GetFileName(indexToAdd);

			/*
			* FILE* betterFile = fopen(fileName.c_str(), "rb");
			if (betterFile != nullptr && betterFile->exists() && betterFile->isFile())
			{
				T read = Read(betterFile);
				if (read != null)
				{
					data.addLast(read);
					if (manager != null)
						manager.RequestCachePresenceAtCurrentResolution(this);
				}
				continue;
			}
			*/
			

			//Second: okay, so we'll have to ask the data provider.
			int betterDim = DataImageDimensions[indexToAdd];
			std::vector<Data> better = Create(betterDim);
			data.addLast(better);
			if (manager != nullptr)
				manager->RequestCachePresenceAtCurrentResolution(this);
			SaveTopResolution();
		}
	}

	std::string GetFileName(int dimIndex)
	{
		if (dimIndex < 0 || dimIndex >= NumDataImageDimensions)
			return "";
		int dim = DataImageDimensions[dimIndex];
		std::string file = "temp";
		//std::string file = parentDir + std::filesystem::path::preferred_separator + dirPrefix;
		//file += dim;
		//file += std::filesystem::path::preferred_separator + folderName + ".rsdat";
		return file;
	}
	bool SaveTopResolution()
	{
		std::vector<Data>& top = data[data.size() - 1];
		std::string fileName = GetFileName(data.size() - 1);
		return SaveImage(top, fileName);
	}

	bool KickTopResolution()
	{
		if (changesRelativeToFile)
			SaveTopResolution();

		std::unique_lock lock(imgLock);
		RemoveBestResolution();
		return true;
	}
private:
	std::string parentDir;
	std::string dirPrefix;
	std::string folderName;
	std::string fullName;
	DataImageManager <Data, ImageType>* manager;
	bool changesRelativeToFile;
	DataProviderBase<ImageType>* source;
	std::shared_mutex imgLock;
	std::vector<std::vector<Data>>data;

};

template <typename Data, typename ImageType>
void
DataImageManager<Data, ImageType>::FullKick(DataImage<Data, ImageType>* requester)
{
	std::unique_lock lock(mutex);
	for (LinkedHashSet<DataImage<Data, ImageType>*> cache : dataImages)
	{
		cache.Remove(requester);
	}
}

template <typename Data, typename ImageType>
void 
DataImageManager<Data, ImageType>::KickingInsert(
	DataImage<Data, ImageType>* add, 
	int index, 
	std::vector<DataImage<Data, ImageType>*> kickRecord)
{
	if (index == -1)
		return;
	dataImages[index].Add(add);
	if(dataImages[index].GetSize() > capacities[index])
	{
		DataImage<Data, ImageType>* kicked = dataImages[index].KickLast();
		kickRecord.push_back(kicked);
		KickingInsert(kicked, index - 1, kickRecord);
	}
}

template <typename Data, typename ImageType>
void 
DataImageManager<Data, ImageType>::RequestCachePresenceAtCurrentResolution(DataImage<Data, ImageType>* requester)
{
	std::shared_lock read(mutex);
	int targetIndex = requester->GetCurrentResolution();
	//We have the image where it should be, we're done
	if (dataImages[targetIndex].Contains(requester))
	{
		return;
	}
	int currentIndex = -1;
	for (int i = NumDataImageDimensions - 1; i >= 0; i--)
	{
		if (dataImages[i].Contains(requester))
		{
			currentIndex = i;
			break;
		}
	}
	//We have the image at better than it needs to be, we're done
	if (currentIndex > targetIndex)
	{
		return;
	}
	read.unlock();

	std::vector<DataImage<Data, ImageType>*> kicked;

	std::unique_lock write(mutex);

	//We need to remove it from its current spot before adding it somewhere else
	if (currentIndex != -1)
	{
		dataImages[currentIndex].Remove(requester);
	}
	KickingInsert(requester, targetIndex, kicked);
	write.unlock();
	for (DataImage<Data, ImageType> di : kicked)
		di.KickTopResolution();
}