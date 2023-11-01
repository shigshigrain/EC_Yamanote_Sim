
#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <thread>



const std::vector<std::vector<double>> shift_wait = {

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅1での調整時間一覧 
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅2での調整時間一覧
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅3
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅4
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅5
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // 駅6以降も記述可能　縦方向の数は駅の数と一致

};


static std::queue<std::vector<double>> wait_que;

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
void make_que_wait() {

	while (!wait_que.empty())wait_que.pop();

	wait_que.emplace(shift_wait.at(0));

	for (size_t i = 1; i < shift_wait.size(); i++) {
		
		while (wait_que.front().size() == i) {

			for (size_t j = 0; j < shift_wait.at(i).size(); j++) {
				std::vector<double> temp = wait_que.front();
				temp.push_back(shift_wait.at(i).at(j));
				wait_que.push(temp);
			}

			wait_que.pop();

		}

	}

}


void simulate(std::vector<double> wait_list) {

	// シミュレーションの中身
	// ここで、駅での調整時間タイミングは以下で取得できる

	// main関数に書いてあったやつをここら辺に上手く貼り付ける。分からなかったら相談

	std::vector<double> wait_time = wait_que.front(); // 先頭のデータを読み取る。次のループでは先頭のデータを捨てることで次のデータを読みだせる
	// ここで読みだした各駅ごとの調整時分を上手くシミュレーションに適応させる


	// 結果をフォルダに出力する  ファイル名をユニークに付けておく

	// シミュレーション本体をmain関数に記述せず関数にしておくと、処理を並列化することができるようになる(ここではしていない)

	return;

}


// main関数もどき
void _main() {

	make_que_wait(); // 調整時分の列挙リスト制作

	std::thread pallarel;

	while (!wait_que.empty()) // リストにシミュレーションすべき時分ずらしパターンが残っていれば、ループ続行
	{
		simulate(wait_que.front()); // 中でwait_queの先頭データを使用して一回分のシミュレーションが終了
		
		//std::vector<double> temp = wait_que.pop(); // ☆シミュレーションし終わった時分ずらしパターンは破棄。つぎのパターンが勝手に先頭になる



	}

	return;

}

