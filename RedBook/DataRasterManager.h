#pragma once

#include <unordered_map>
#include <unordered_set>
#include <list>
#include <array>
#include <stdexcept>
#include <memory>
#include "UniqueRequestStack.h"

class WorldMap;

struct DataRasterDefaultPolicy
{
	static bool ThreadForAllocation() { return false; }

	template <typename Handle>
	static void OnAllocation(Handle handle, WorldMap* world) {}

	template <typename Handle>
	static void OnDeallocation(Handle handle, WorldMap* world) {}
};

template <typename DataType, int w, int h, typename ID, typename IDHash, typename Policy = DataRasterDefaultPolicy>
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
			return data[w * y + x];
		}
		inline DataType Get(uint32_t index)
		{
			return data[index];
		}
		inline void Set(DataType d, uint32_t x, uint32_t y)
		{
			data[w * y + x] = d;
		}
		inline void Set(DataType d, uint32_t index)
		{
			data[index] = d;
		}
		inline void Clear()
		{
			std::fill(std::begin(data), std::end(data), DataType{});
		}
		DataType data[w * h] = {};
	};

//------------------------------------------------------------------------------

	struct DataRasterSlot
	{
		DataRasterSlot(ID assigned) : assignedID(assigned)
		{
			dataReady = false;
		}
		~DataRasterSlot()
		{

		}
		void Invalidate()
		{
			dataReady = false;
		}
		void MarkDataReady()
		{
			dataReady = true;
		}
		void Reassign(ID newID)
		{
			dataReady = false;
			assignedID = newID;
			data.Clear();
		}
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
		bool OkayToUse() const
		{
			return Valid() && Ready();
		}
		bool Valid() const
		{
			return data != nullptr && myID == data->assignedID;
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
		inline DataType* GetRawData()
		{
			if (!Valid())
				throw std::runtime_error("Dereferencing an outdated Data Raster Handle");
			return data->data.data;
		}

	protected:
		void UpdateSlot(SlotPtr newSlot)
		{
			data = newSlot;
		}
	public:
		ID myID;
	private:
		friend class DataRasterManager;
		LimitedDataRasterHandle(SlotPtr slotToAccept, ID myID) : data(slotToAccept), myID(myID)
		{
		}
		SlotPtr data;
	};
	
//------------------------------------------------------------------------------
// 	   MEMBER VARIABLES
//------------------------------------------------------------------------------

private:
	int fullCapacity;

	std::unordered_map<SlotPtr, typename std::list<SlotPtr>::iterator> locLookup;
	std::list<SlotPtr> cache;
	std::unordered_map<ID, SlotPtr, IDHash> idLookup;
	std::unordered_set<SlotPtr> pinned;

	RequestStack<ID, IDHash>* threads;

	std::mutex manager_mutex;

	WorldMap* myWorld;

private:

//------------------------------------------------------------------------------
// 	   METHODS NEEDED FOR PUBLIC HANDLE CLASS
//------------------------------------------------------------------------------

	SlotPtr GetSlotForID(const ID& id)
	{
		std::unique_lock<std::mutex> lock(manager_mutex);
		auto it = idLookup.find(id);
		if (it != idLookup.end())
		{
			// We already have this raster in cache
			SlotPtr slot = it->second;
			if (pinned.find(slot) != pinned.end())
			{
				//we have it, it's pinned, we're fine
				return slot;
			}

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
			return slot;
		}
		auto allocationDecision = [this, id]()
		{
			std::unique_lock<std::mutex> innerLock(manager_mutex);
			SlotPtr slotToUse = nullptr;
			if (cache.size() + pinned.size() < fullCapacity)
			{
				SlotPtr newSlot = SlotPtr(new DataRasterSlot(id));
				pinned.insert(newSlot);
				idLookup[id] = newSlot;
				slotToUse = newSlot;
				innerLock.unlock();
			}
			else
			{
				if (cache.empty())
					return slotToUse; //returns nullptr

				DataRasterSlot* last = cache.back();
				ID prevID = last->assignedID;
				cache.pop_back();
				locLookup.erase(last);
				pinned.insert(last);

				innerLock.unlock();
				Policy::OnDeallocation(LimitedDataRasterHandle(last, prevID), myWorld);
				innerLock.lock();

				last->Invalidate();
				idLookup.erase(prevID);
				last->Reassign(id);
				idLookup[id] = last;
				slotToUse = last;
				innerLock.unlock();
			}

			Policy::OnAllocation(LimitedDataRasterHandle(slotToUse, id), myWorld);

			innerLock.lock();
			pinned.erase(slotToUse);
			cache.push_front(slotToUse);
			locLookup[slotToUse] = cache.begin();
			slotToUse->MarkDataReady();

			return slotToUse;
		};

		if (threads != nullptr && Policy::ThreadForAllocation())
		{
			threads->RequestTask(id, allocationDecision);
			return nullptr;
		}
		else
			return allocationDecision();
	}

	//------------------------------------------------------------------------------


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
	DataRasterManager(int capacity, WorldMap* myWorld, RequestStack<ID, IDHash>* pool = nullptr) : myWorld(myWorld), threads(pool)
	{
		if (capacity <= 0)
			throw std::runtime_error("DataRasterManager must have a positive capacity");
		fullCapacity = capacity;
	}

	void RefreshRaster(const ID& id)
	{
		GetSlotForID(id);
	}

	DataRasterHandle GetRaster(const ID& id)
	{
		return DataRasterHandle(GetSlotForID(id), id, this);
	}

	int GetMaxCapacity() const { return fullCapacity; }
	int GetNumUsed() const { return (int)(cache.size() + pinned.size()); }
};