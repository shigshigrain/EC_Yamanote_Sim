#pragma once

#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <mutex>
#include <chrono>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <filesystem>
#include <source_location>

class CsvWriter
{

public:

	std::ofstream fp7;
	std::string FileNameH;
	std::queue<std::string> TaskQue; // 非同期タスクとしてアクセス

public:

	CsvWriter();
	bool Init();
	bool Init(std::string header);
	bool StandUp();
	bool StandUp(std::string header);
	bool PushTask(std::string str, std::mutex& _mtx);
	bool Writing(std::mutex& _mtx);

	~CsvWriter();

};

std::string getTimeStamp();