#include "CsvWriter.hpp"

CsvWriter::CsvWriter()
{

	Init();

}

bool CsvWriter::Init()
{

	FileNameH = "";

	return false;
}

bool CsvWriter::Init(std::string header)
{
	FileNameH = header;

	return true;
}

bool CsvWriter::StandUp()
{

	std::string file_name("Result\\2列車_消費電力_山手線_5駅5パターン_" + getTimeStamp() + ".csv");

	try
	{
		fp7 = std::ofstream(file_name);
		fp7.exceptions(std::ios_base::failbit);
	}
	catch (const std::exception& e) {
		const std::source_location location = std::source_location::current();
		std::cerr << "ファイルを開けませんでした。 #CsvWriter.cpp <L" << location.line() << "> #" << e.what() << std::endl;
		return false;
	}

	return true;
}

bool CsvWriter::StandUp(std::string header)
{
	FileNameH = header;

	std::string file_name("Result\\" + FileNameH +"内1外1列車_消費電力_5駅4パターン_" + getTimeStamp() + ".csv");

	try
	{
		fp7 = std::ofstream(file_name);
		fp7.exceptions(std::ios_base::failbit);
	}
	catch (const std::exception& e) {
		const std::source_location location = std::source_location::current();
		std::cerr << "ファイルを開けませんでした。 #CsvWriter.cpp <L" << location.line() << "> #" << e.what() << std::endl;
		return false;
	}

	return true;
}

bool CsvWriter::PushTask(std::string str, std::mutex& _mtx)
{

	std::lock_guard<std::mutex> _lock(_mtx);

	TaskQue.push(str);

	return true;
}

bool CsvWriter::Writing(std::mutex& _mtx)
{
	std::lock_guard<std::mutex> _lock(_mtx);

	while (!TaskQue.empty()) {

		fp7 << TaskQue.front() << "\n";
		TaskQue.pop();

	}

	return true;
}

CsvWriter::~CsvWriter()
{

	//fp7.close();


}

std::string getTimeStamp()
{
	using namespace std;
	#define tos to_string

	time_t start; struct tm today;
	time(&start);
	errno_t error = localtime_s(&today, &start);

	return string(tos(today.tm_year) + tos(today.tm_mon + 1) + tos(today.tm_mday) + tos(today.tm_hour) + tos(today.tm_min) + tos(today.tm_sec));
}