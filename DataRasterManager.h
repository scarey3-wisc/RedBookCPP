#pragma once

#include <unordered_map>
#include <unordered_set>
#include <list>
#include <array>
#include <stdexcept>
#include <memory>

class WorldMap;

struct DataRasterDefaultPolicy
{
	static bool ThreadOnAllocation() { return false; }
	static bool ThreadOnDeallocation() { return false; }

	template <typename Handle>
	static void OnAllocation(Handle handle, WorldMap* world) {}

	template <typename Handle>
	static void OnDeallocation(Handle handle, WorldMap* world) {}
};

template <typename DataType, int dim, typename ID, typename IDHash, typename Policy = DataRasterDefaultPolicy>
class DataRasterManager
{
//------------------------------------------------------------------------------
// 	   SIMPLE INNER CLASSES
//------------------------------------------------------------------------------

private:

	struct alignas(64) Raster
	{
		inline DataType Get(uint32_t x, uint32_t y)
		{
			return data[dim * y + x];
		}
		inline DataType Get(uint32_t index)
		{
			return data[index];
		}
		inline void Set(DataType d, uint32_t x, uint32_t y)
		{
			data[dim * y + x] = d;
		}
		inline void Set(DataType d, uint32_t index)
		{
			data[index] = d;
		}
		inline void Clear()
		{
			std::fill(std::begin(data), std::end(data), DataType{});
		}
		DataType data[dim * dim] = {};
	};

//------------------------------------------------------------------------------

	struct DataRasterSlot
	{
		DataRasterSlot(ID assigned) : assignedID(assigned)
		{
			useID = 0;
			dataReady = false;
		}
		~DataRasterSlot()
		{

		}
		void MarkDataReady()
		{
			dataReady = true;
		}
		void Reassign(ID newID)
		{
			useID++;
			dataReady = false;
			assignedID = newID;
			data.Clear();
		}
		int useID;
		ID assignedID;
		bool dataReady;
		Raster data;
	};

//------------------------------------------------------------------------------

	using SlotPtr = DataRasterSlot*;

//------------------------------------------------------------------------------

public:
	class LimitedDataRasterHandle
	{
	public:
		LimitedDataRasterHandle(const LimitedDataRasterHandle& other) = default;
		LimitedDataRasterHandle& operator=(const LimitedDataRasterHandle& other) = default;
		bool Valid() const
		{
			return useID == data->useID;
		}
		bool Ready() const
		{
			if (!Valid())
				throw std::runtime_error("Checking readiness of an outdated Data Raster Handle");
			return data->dataReady;
		}
		inline DataType Get(uint32_t x, uint32_t y)
		{
			if (!Valid())
				throw std::runtime_error("Dereferencing an outdated Data Raster Handle");
			return data->data.Get(x, y);
		}
		inline void Set(DataType d, uint32_t x, uint32_t y)
		{
			if (!Valid())
				throw std::runtime_error("Dereferencing an outdated Data Raster Handle");
			data->data.Set(d, x, y);
		}
		inline DataType Get(uint32_t index)
		{
			if (!Valid())
				throw std::runtime_error("Dereferencing an outdated Data Raster Handle");
			return data->data.Get(index);
		}
		inline void Set(DataType d, uint32_t index)
		{
			if (!Valid())
				throw std::runtime_error("Dereferencing an outdated Data Raster Handle");
			data->data.Set(d, index);
		}

	protected:
		void UpdateSlot(SlotPtr newSlot)
		{
			data = newSlot;
			useID = data->useID;
		}
	public:
		ID myID;
	private:
		friend class DataRasterManager;
		LimitedDataRasterHandle(SlotPtr slotToAccept, ID myID) : data(slotToAccept), myID(myID)
		{
			useID = data->useID;
		}
		int useID;
		SlotPtr data;
	};
	
//------------------------------------------------------------------------------
// 	   MEMBER VARIABLES
//------------------------------------------------------------------------------

private:
	int fullCapacity;
	int currentCacheCapacity;

	std::unordered_map<SlotPtr, typename std::list<SlotPtr>::iterator> locLookup;
	std::list<SlotPtr> cache;
	std::unordered_map<ID, SlotPtr, IDHash> idLookup;
	std::unordered_set<SlotPtr> pinned;

	WorldMap* myWorld;

private:

//------------------------------------------------------------------------------
// 	   METHODS NEEDED FOR PUBLIC HANDLE CLASS
//------------------------------------------------------------------------------

	SlotPtr UnpinSlotForID(const ID& id)
	{
		auto it = idLookup.find(id);
		if (it != idLookup.end())
		{
			// We do have this raster in cache
			SlotPtr slot = it->second;
			auto pinIt = pinned.find(slot);
			if (pinIt != pinned.end())
			{
				pinned.erase(pinIt);
				// Move to front of cache
				cache.push_front(slot);
				locLookup[slot] = cache.begin();
				currentCacheCapacity++;
				return slot;
			}
			else
			{
				// It wasn't pinned
				return nullptr;
			}
		}
		else
		{
			// We don't have this raster at all
			return nullptr;
		}
	}

	SlotPtr GetPinnedSlotForID(const ID& id)
	{
		auto it = idLookup.find(id);
		if (it != idLookup.end())
		{
			// We already have this raster in cache
			SlotPtr slot = it->second;
			if (pinned.find(slot) == pinned.end())
			{
				auto locIt = locLookup.find(slot);
				if (locIt != locLookup.end())
				{
					// Remove from cache, add to pinned list
					currentCacheCapacity--;
					cache.erase(locIt->second);
					locLookup.erase(locIt->first);
					pinned.insert(slot);
				}
				else
				{
					throw new std::runtime_error("DataRasterManager internal error: we claim to be managing a slot but have no record of its location");
				}
			}
			return slot;
		}

		if (cache.size() >= currentCacheCapacity)
		{
			// We need to kick something out of the cache
			auto last = cache.back();
			//TODO: we need to give last an opportunity to save itself if it's dirty
			Policy::OnDeallocation(LimitedDataRasterHandle(last, last->assignedID), myWorld);

			//Now we can remove last from the cache.
			currentCacheCapacity--;
			cache.pop_back();
			idLookup.erase(last->assignedID);
			locLookup.erase(last);

			//TODO: we need to populate last with the new data for ID here
			last->Reassign(id);
			Policy::OnAllocation(LimitedDataRasterHandle(last, id), myWorld);
			last->MarkDataReady();
			//And now we reassign it to the new ID and add it to the pinned area
			idLookup[id] = last;
			pinned.insert(last);
			return last;
		}

		else
		{
			// We have space to make a new slot
			SlotPtr newSlot = SlotPtr(new DataRasterSlot(id));
			Policy::OnAllocation(LimitedDataRasterHandle(newSlot, id), myWorld);
			newSlot->MarkDataReady();
			currentCacheCapacity--;
			idLookup[id] = newSlot;
			pinned.insert(newSlot);
			return newSlot;
		}
	}

	SlotPtr GetSlotForID(const ID& id)
	{
		auto it = idLookup.find(id);
		if (it != idLookup.end())
		{
			// We already have this raster in cache
			SlotPtr slot = it->second;
			if (pinned.find(slot) == pinned.end())
			{
				auto locIt = locLookup.find(slot);
				if (locIt != locLookup.end())
				{
					// Move to front of cache
					cache.splice(cache.begin(), cache, locIt->second);
				}
				else
				{
					throw new std::runtime_error("DataRasterManager internal error: we claim to be managing a slot but have no record of its location");
				}
			}
			return slot;
		}

		if (cache.size() >= currentCacheCapacity)
		{
			// We need to kick something out of the cache
			auto last = cache.back();
			//TODO: we need to give last an opportunity to save itself if it's dirty
			Policy::OnDeallocation(LimitedDataRasterHandle(last, last->assignedID), myWorld);

			//Now we can remove last from the cache.
			cache.pop_back();
			idLookup.erase(last->assignedID);
			locLookup.erase(last);

			//TODO: we need to populate last with the new data for ID here
			last->Reassign(id);
			Policy::OnAllocation(LimitedDataRasterHandle(last, id), myWorld);
			last->MarkDataReady();
			//And now we reassign it to the new ID and add it to the front of the cache
			idLookup[id] = last;
			cache.push_front(last);
			locLookup[last] = cache.begin();
			return last;
		}

		else
		{
			// We have space to make a new slot
			SlotPtr newSlot = SlotPtr(new DataRasterSlot(id));
			Policy::OnAllocation(LimitedDataRasterHandle(newSlot, id), myWorld);
			newSlot->MarkDataReady();
			idLookup[id] = newSlot;
			cache.push_front(newSlot);
			locLookup[newSlot] = cache.begin();
			return newSlot;
		}
	}

public:

//------------------------------------------------------------------------------
// 	   PUBLIC HANDLE CLASS
//------------------------------------------------------------------------------

	class DataRasterHandle : public LimitedDataRasterHandle
	{
	public:
		DataRasterHandle(const DataRasterHandle& other) = default;
		DataRasterHandle& operator=(const DataRasterHandle& other) = default;
		void Refresh()
		{
			LimitedDataRasterHandle::UpdateSlot(source->GetSlotForID(LimitedDataRasterHandle::myID));
		}
		void Pin()
		{
			LimitedDataRasterHandle::UpdateSlot(source->GetPinnedSlotForID(LimitedDataRasterHandle::myID));
		}
		void Unpin()
		{
			SlotPtr newData = source->UnpinSlotForID(LimitedDataRasterHandle::myID);
			if (newData == nullptr)
				throw std::runtime_error("Attempt to unpin a Data Raster Handle that wasn't pinned");
			LimitedDataRasterHandle::UpdateSlot(newData);
		}

	private:
		DataRasterHandle(SlotPtr slotToAccept, ID myID, DataRasterManager* source) : 
			LimitedDataRasterHandle(slotToAccept, myID), source(source)
		{
			
		}
		friend class DataRasterManager;
		DataRasterManager* source;
	};

//------------------------------------------------------------------------------
// 	   FULLY DEFINED METHODS
//------------------------------------------------------------------------------

public:
	DataRasterManager(int capacity, WorldMap* myWorld) : myWorld(myWorld)
	{
		if (capacity <= 0)
			throw std::runtime_error("DataRasterManager must have a positive capacity");
		fullCapacity = capacity;
		currentCacheCapacity = capacity;
	}

	DataRasterHandle GetRaster(const ID& id)
	{
		return DataRasterHandle(GetSlotForID(id), id, this);
	}
};