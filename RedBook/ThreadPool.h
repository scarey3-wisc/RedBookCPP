#pragma once

#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>

class ThreadPool {
public:
	// // Constructor to creates a thread pool with given
	// number of threads
	ThreadPool(size_t num_threads = 1)
	{

		// Creating worker threads
		for (size_t i = 0; i < num_threads; ++i) {
			threads.emplace_back([this] {
				while (true) {
					std::function<void()> task;
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
							return !tasks.empty() || stop;
							});

						// exit the thread in case the pool
						// is stopped and there are no tasks
						if (stop && tasks.empty()) {
							return;
						}

						// Get the next task from the queue
						task = std::move(tasks.front());
						tasks.pop();
					}

					task();
				}
				});
		}
	}

	// Destructor to stop the thread pool
	~ThreadPool()
	{
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

	// Enqueue task for execution by the thread pool
	void enqueue(std::function<void()> task)
	{
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			tasks.emplace(move(task));
		}
		cv.notify_one();
	}
	void discardAllQueuedTasks()
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		while (!tasks.empty())
			tasks.pop();
	}
	size_t NumTasksQueued()
	{
		return tasks.size();
	}

private:
	// Vector to store worker threads
	std::vector<std::thread> threads;

	// Queue of tasks
	std::queue<std::function<void()>> tasks;

	// Mutex to synchronize access to shared data
	std::mutex queue_mutex;

	// Condition variable to signal changes in the state of
	// the tasks queue
	std::condition_variable cv;

	// Flag to indicate whether the thread pool should stop
	// or not
	bool stop = false;
};