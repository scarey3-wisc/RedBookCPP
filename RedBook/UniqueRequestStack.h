#pragma once

#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <map>

template <typename ID, typename IDHash>
struct RequestStack {
private:
	int CheckPriority(ID id)
	{
		return id.GetPriority();
	}
public:
	// // Constructor to creates a thread pool with given
	// number of threads
	RequestStack(size_t num_threads = 1)
	{

		// Creating worker threads
		for (size_t i = 0; i < num_threads; ++i) {
			threads.emplace_back([this] {
				while (true) {
					std::function<void()> task;
					ID taskID;
					// The reason for putting the below code
					// here is to unlock the queue before
					// executing the task so that other
					// threads can perform enqueue tasks
					{
						// Locking the queue so that data
						// can be shared safely
						std::unique_lock<std::mutex> lock(queue_mutex);

						// Waiting until there is a task to
						// execute or the pool is stopped
						cv.wait(lock, [this] {
							if (stop)
								return true;
							bool allEmpty = true;
							for (auto& stack : taskStack)
								if (stack.second.size() > 0)
								{
									allEmpty = false;
									break;
								}
							return !allEmpty || stop;
							});

						// exit the thread in case the pool
						// is stopped and there are no tasks
						bool allEmpty = true;
						for (auto& stack : taskStack)
							if (stack.second.size() > 0)
							{
								allEmpty = false;
								break;
							}
						if (stop && allEmpty) 
						{
							return;
						}

						// Get the next task from the queue

						std::pair< ID, std::function<void()>> request = std::move(taskStack.begin()->second.front());
						task = std::move(request.second);
						taskID = std::move(request.first);
						taskStack.begin()->second.pop_front();
						requestLocation.erase(request.first);
						if(taskStack.begin()->second.size() == 0)
							taskStack.erase(taskStack.begin());
					}

					task();
					{
						std::unique_lock<std::mutex> lock(queue_mutex);
						requestLookup.erase(taskID);
					}
				}
				});
		}
	}

	void TerminateThreads()
	{
		if(stop)
			return;
		{
			// Lock the queue to update the stop flag safely
			std::unique_lock<std::mutex> lock(queue_mutex);
			stop = true;
		}

		// Notify all threads
		cv.notify_all();

		// Joining all worker threads to ensure they have
		// completed their tasks
		for (auto& thread : threads) {
			thread.join();
		}
	}

	// Destructor to stop the thread pool
	~RequestStack()
	{
		TerminateThreads();
	}

	// Enqueue task for execution by the thread pool
	void RequestTask(ID id, std::function<void()> task)
	{
		int prio = CheckPriority(id);
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			if (requestLookup.find(id) != requestLookup.end())
			{
				auto locIt = requestLocation.find(id);
				if (locIt != requestLocation.end())
					taskStack[prio].splice(taskStack[prio].begin(), taskStack[prio], locIt->second);
				return;
			}
			requestLookup.insert(id);
			taskStack[prio].push_front({id, task});
			requestLocation[id] = taskStack[prio].begin();
		}
		cv.notify_one();
	}
	void DiscardAllWaitingTasks()
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		//this clearing pattern might seem strange when we could just clear taskStack and
		//requestLocation, but requestLookup can contain IDs that aren't in either (because
		//a thread has started working on them), so we need to be careful here to properly
		//clear them
		for (auto task : requestLocation)
		{
			ID taskID = task.first;
			int taskPriority = CheckPriority(taskID);
			taskStack[taskPriority].erase(task.second);
			if (taskStack[taskPriority].size() == 0)
				taskStack.erase(taskPriority);
			requestLookup.erase(task.first);
		}
		requestLocation.clear();
	}
	size_t NumActiveRequests()
	{
		return requestLookup.size();
	}

private:
	// Vector to store worker threads
	std::vector<std::thread> threads;

	// Queue of tasks
	std::map<int, std::list<std::pair<ID, std::function<void()>>>> taskStack;
	std::unordered_set<ID, IDHash> requestLookup;
	std::unordered_map<ID, typename std::list<std::pair<ID, std::function<void()>>>::iterator, IDHash> requestLocation;

	// Mutex to synchronize access to shared data
	std::mutex queue_mutex;

	// Condition variable to signal changes in the state of
	// the tasks queue
	std::condition_variable cv;

	// Flag to indicate whether the thread pool should stop
	// or not
	bool stop = false;
};