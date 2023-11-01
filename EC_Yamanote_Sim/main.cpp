
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


const std::vector<std::vector<double>> shift_wait = {

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅1での調整時間一覧 
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅2での調整時間一覧
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅3
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅4
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅5
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅6以降も記述可能　縦方向の数は駅の数と一致

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
// が列挙されている状況になっているはず

// 調整時分パターンの列挙リスト作成
void make_que_wait(std::queue<std::vector<double>>& wq) {

	while (!wq.empty())wq.pop();

	for (auto&& _sw : shift_wait.at(0)) {
		std::vector<double> _tmp(1, _sw);
		wq.push(_tmp);
	}
	
	for (size_t i = 1; i < shift_wait.size(); i++) {

		while (wq.front().size() == i) {

			for (size_t j = 0; j < shift_wait.at(i).size(); j++) {
				std::vector<double> temp = wq.front();
				temp.push_back(shift_wait.at(i).at(j));
				wq.push(temp);
			}

			wq.pop();

		}

	}

}

constexpr size_t NumThread = 3; // 

void simulate() {
	
	std::queue<std::vector<double>> wait_que;

	make_que_wait(wait_que); // 調整時分の列挙リスト制作
	size_t counter = 0;

	std::mutex mtx;

	CsvWriter CW;
	CW.StandUp("Tester_");
	CW.fp7 << "目白発車時刻変化分[sec],高田馬場[sec],新大久保[sec],代々木[sec],原宿[sec],変電所出力総電力量[kWh]\n";

	std::vector<std::thread> TaskPool(0);
	TaskPool.reserve(NumThread);
	std::vector<SimulateEngine> se(NumThread);

	std::cout << "Starting Simulation \n";

	while (!wait_que.empty()) // リストにシミュレーションすべき時分ずらしパターンが残っていれば、ループ続行
	{
		
		//se.clear();
		//se = std::vector<SimulateEngine>(NumThread);

		for (size_t ii = 0; ii < NumThread; ii++) {
			
			if (wait_que.empty())break;
			
			std::vector<double> wf = wait_que.front();
			//TaskPool.push_back(std::thread(mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx))));
			//TaskPool.emplace_back(std::thread([wait_que, &CW, &mtx]() {mySimulate(wait_que.front(), std::ref(CW), std::ref(mtx)); }));
			TaskPool.emplace_back(std::thread([ii, &se, wf, &CW, &mtx]() {ExeSimulate(std::ref(se.at(ii)), wf, std::ref(CW), std::ref(mtx)); }));
			wait_que.pop();

		}

		for (size_t i = 0; i < TaskPool.size(); i++) {
			TaskPool.at(i).join(); // 分解したスレッドの処理すべてをここでいったん待つ
		}

		TaskPool.clear();

		//se.clear();

		counter += NumThread;

		std::cout << "Finish " << counter << " Thread\n";

		if (counter >= 50Ui64)break;

	}

	CW.Writing(mtx);
	std::cout << "Finish Writing \n";

	CW.~CsvWriter();

	std::string waiter;
	std::cin >> waiter;

	return;

}


// main関数もどき
int main(void) {
	
	simulate();

	std::cout << "\n Finish Simulating\n Type to End\n";

	std::string waiter;
	std::cin >> waiter;

	return 0;

}
