
#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <thread>
#include <future>
#include <memory>

//#include "mySimulate.h"
//#include "CsvWriter.hpp"

#include "SimulateEngine.hpp"


const std::vector<std::vector<double>> ShiftWait5 = {

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅1での調整時間一覧^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅2^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅3^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅4^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅5^
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅6以降も記述可能　縦方向の数は駅の数と一致^
};

const std::vector<std::vector<double>> ShiftWaitIN = {

	{-5.0, -4.0, 0.0}, // 目白^
	{-5.0, -4.0, 0.0}, // 高田馬場^
	{-5.0, -4.0, 0.0}, // 新大久保^
	{-5.0, -4.0, 0.0}, // 代々木^
	{-5.0, -4.0, 0.0}, // 原宿^
};

const std::vector<std::vector<double>> ShiftWaitOUT = {

	{-5.0, -4.0, 0.0}, // 原宿^
	{-5.0, -4.0, 0.0}, // 代々木^
	{-5.0, -4.0, 0.0}, // 新大久保^
	{-5.0, -4.0, 0.0}, // 高田馬場^
	{-5.0, -4.0, 0.0}, // 目白^
};

// 中身は
/*
{ 駅1の調整時分, 駅2の調整時分, 駅3・・・, 駅4・・・, 駅5・・・} // パターン1
{ 駅1の調整時分, 駅2の調整時分, 駅3・・・, 駅4・・・, 駅5・・・} // パターン2
						・
						・
						・
// 具体的に
{ -5.0, -5.0, -5.0, -5.0, -5.0 } // パターン1
{ -5.0, -5.0, -5.0, -5.0, -2.5 } // パターン2
{ -5.0, -5.0, -5.0, -5.0, 0.0 } // パターン3
{ -5.0, -5.0, -5.0, -5.0, 2.5 } // パターン4
{ -5.0, -5.0, -5.0, -5.0, 5.0 } // パターン5
{ -5.0, -5.0, -5.0, -2.5, -5.0 } // パターン6
{ -5.0, -5.0, -5.0, -2.5, -2.5 } // パターン7
{ -5.0, -5.0, -5.0, -2.5, 0.0 } // パターン8
						・
						・
						・
{ 5.0, 5.0, 5.0, 5.0, 2.5 } // パターン3124
{ 5.0, 5.0, 5.0, 5.0, 5.0 } // パターン3125
*/
// が列挙されている状況になっているはず^

// 調整時分パターンの列挙リスト作成^
static void make_que_wait(std::deque<std::vector<double>>& wq, const std::vector<std::vector<double>>& SeedF) {

	while (!wq.empty())wq.pop_back();

	for (auto&& _sw : SeedF.at(0)) {
		std::vector<double> _tmp(1, _sw);
		wq.push_back(_tmp);
	}
	
	for (size_t i = 1; i < SeedF.size(); i++) {

		while (wq.front().size() == i) {

			for (size_t j = 0; j < SeedF.at(i).size(); j++) {
				std::vector<double> temp = wq.front();
				temp.push_back(SeedF.at(i).at(j));
				wq.push_back(temp);
			}

			wq.pop_front();

		}

	}

}

constexpr size_t NumThread = 3; // 

static void simulate() {
	
	std::deque<std::vector<double>> wait_IN;
	std::deque<std::vector<double>> wait_OUT;

	make_que_wait(wait_IN, ShiftWaitIN); // 調整時分の列挙リスト制作^
	make_que_wait(wait_OUT, ShiftWaitOUT); // 調整時分の列挙リスト制作^
	size_t counter = 0;

	std::mutex mtx;

	CsvWriter CW;
	CW.StandUp("Tester_");
	CW.fp7 << "内:目白停車時刻変化[s],内:高田馬場[s],内:新大久保[s],内:代々木[s],内:原宿[s],"
		   << "外:原宿停車時刻変化[s],外:代々木[s],外:新大久保[s],外:高田馬場[s],外:目白[s],"
		   << "変電所出力総電力量[kWh]\n";

	std::vector<std::thread> TaskPool(0);
	TaskPool.reserve(NumThread);
	std::vector<SimulateEngine> se(NumThread);

	std::cout << "Starting Simulation \n";

	// 内回り・外回りそれぞれのループで実装^
	for (size_t i = 0; i < wait_IN.size(); i++) {
		size_t j = 0;
		while (j < wait_OUT.size())
		{

			for (size_t ii = 0; ii < NumThread; ii++) {

				if (j >= wait_OUT.size())break;
				
				std::vector<double> wf_IN = wait_IN.at(i);
				std::vector<double> wf_OUT = wait_OUT.at(j);
				
				//TaskPool.push_back(std::thread(mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx))));
				//TaskPool.emplace_back(std::thread([wait_que, &CW, &mtx]() {mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx)); }));
				TaskPool.emplace_back(std::thread([ii, &se, wf_IN, wf_OUT, &CW, &mtx]() {ExeSimulate(std::ref(se.at(ii)), wf_IN, wf_OUT, std::ref(CW), std::ref(mtx)); }));
				j++;
			}

			for (auto&& tp : TaskPool) {
				tp.join(); // 分解したスレッドの処理すべてをここでいったん待つ
			}

			TaskPool.clear();

			//se.clear();

			counter += NumThread;

			std::cout << "Finish " << counter << " Thread\n";

			//if (j >= 2Ui64)break;

			
		}


		CW.Writing(mtx);
		if (i >= 2Ui64)break;

	}


	//while (!wait_que.empty()) // リストにシミュレーションすべき時分ずらしパターンが残っていれば、ループ続行
	//{
	//	
	//	//se.clear();
	//	//se = std::vector<SimulateEngine>(NumThread);
	//	for (size_t ii = 0; ii < NumThread; ii++) {
	//		
	//		if (wait_que.empty())break;
	//		
	//		std::vector<double> wf = wait_que.front();
	//		//TaskPool.push_back(std::thread(mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx))));
	//		//TaskPool.emplace_back(std::thread([wait_que, &CW, &mtx]() {mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx)); }));
	//		TaskPool.emplace_back(std::thread([ii, &se, wf, &CW, &mtx]() {ExeSimulate(std::ref(se.at(ii)), wf, std::ref(CW), std::ref(mtx)); }));
	//		wait_que.pop();
	//	}
	//	for (size_t i = 0; i < TaskPool.size(); i++) {
	//		TaskPool.at(i).join(); // 分解したスレッドの処理すべてをここでいったん待つ
	//	}
	//	TaskPool.clear();
	//	//se.clear();
	//	counter += NumThread;
	//	std::cout << "Finish " << counter << " Thread\n";
	//	if (counter >= 50Ui64)break;
	//}


	//CW.Writing(mtx);
	std::cout << "Finish Writing \n";

	CW.fp7.close();

	//std::string waiter;
	//std::cin >> waiter;

	se.clear();

	return;

}


int main(void) {
	
	simulate();

	std::cout << "\n Finish Simulating\n Type to End\n";

	//std::string waiter;
	//std::cin >> waiter;

	return 0;

}
