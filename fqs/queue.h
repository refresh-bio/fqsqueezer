#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

// Generic multithreading queues

#include <queue>
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>

using namespace std;


// *****************************************************************************************
//
class CBarrier
{
public:
	CBarrier(const CBarrier&) = delete;
	CBarrier& operator=(const CBarrier&) = delete;
	explicit CBarrier(unsigned int count) :
		m_count(count), m_generation(0),
		m_count_reset_value(count)
	{
	}
	void count_down_and_wait()
	{
		std::unique_lock< std::mutex > lock(m_mutex);
		unsigned int gen = m_generation;
		if (--m_count == 0)
		{
			m_generation++;
			m_count = m_count_reset_value;
			m_cond.notify_all();
			return;
		}
		while (gen == m_generation)
			m_cond.wait(lock);
	}
private:
	std::mutex m_mutex;
	std::condition_variable m_cond;
	unsigned int m_count;
	unsigned int m_generation;
	unsigned int m_count_reset_value;
};

// *****************************************************************************************
//
class CSemaphore {
protected:
	int counter;
	int generation;
	std::mutex mutex;
	std::condition_variable cv;

public:
	// *****************************************************************************************
	//
	CSemaphore(int _counter = 0) : counter(_counter), generation(0) {}

	// *****************************************************************************************
	//
	void Inc(int new_generation = 0) {
		std::unique_lock<std::mutex> lk(mutex);
		if(generation == new_generation)
			++counter;
		else
		{
			generation = new_generation;
			counter = 1;
		}
	}

	// *****************************************************************************************
	//
	void IncNum(int num, int new_generation = 0) {
		std::unique_lock<std::mutex> lk(mutex);
		if(generation == new_generation)
			counter += num;
		else
		{
			generation = new_generation;
			counter = num;
		}
	}

	// *****************************************************************************************
	//
	void Dec(int dec_generation = 0) {
		std::unique_lock<std::mutex> lk(mutex);
		
		if(generation == dec_generation)
			--counter;

		if (counter == 0)
			cv.notify_one();
	}

	// *****************************************************************************************
	//
	void Dec_notify_all(int dec_generation = 0) {
		std::unique_lock<std::mutex> lk(mutex);
		if(generation == dec_generation)
			--counter;

		if(counter == 0)
			cv.notify_all();
	}

	// *****************************************************************************************
	//
	void WaitForZero(int wait_generation = 0) {
		std::unique_lock<std::mutex> lk(mutex);
		cv.wait(lk, [this, &wait_generation] {
			return (counter == 0) && (generation == wait_generation);});
	}
};

// ************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CRegisteringQueue
{
	typedef queue<T, deque<T>> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	// *****************************************************************************************
	//
	CRegisteringQueue(int _n_producers)
	{
		Restart(_n_producers);
	};

	// *****************************************************************************************
	//
	~CRegisteringQueue()
	{};

	// *****************************************************************************************
	//
	void Restart(int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers  = _n_producers;
		n_elements = 0;
	}

	// *****************************************************************************************
	//
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *****************************************************************************************
	//
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *****************************************************************************************
	//
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if(!n_producers)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void Push(T data)
	{
		unique_lock<mutex> lck(mtx);
		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		if(was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void PushRange(vector<T> &data)
	{
		unique_lock<mutex> lck(mtx);
		bool was_empty = n_elements == 0;
		
		for(auto &x : data)
			q.push(x);
		n_elements += data.size();

		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	bool Pop(T &data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_producers;}); 

		if(n_elements == 0)
			return false;

		data = q.front();
		q.pop();
		--n_elements;
		if(n_elements == 0)
			cv_queue_empty.notify_all();

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
		return n_elements;
	}
};

// EOF
