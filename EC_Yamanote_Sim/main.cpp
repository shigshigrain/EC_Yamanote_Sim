
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

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w1�ł̒������Ԉꗗ 
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w2�ł̒������Ԉꗗ
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w3
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w4
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w5
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w6�ȍ~���L�q�\�@�c�����̐��͉w�̐��ƈ�v

};

// ���g��
/*
{ �w1�̒�������, �w2�̒�������, �w3�E�E�E, �w4�E�E�E, �w5�E�E�E} // �p�^�[��1
{ �w1�̒�������, �w2�̒�������, �w3�E�E�E, �w4�E�E�E, �w5�E�E�E} // �p�^�[��2
						�E
						�E
						�E
// ��̓I��
{ -5.0, -5.0, -5.0, -5.0, -5.0 } // �p�^�[��1
{ -5.0, -5.0, -5.0, -5.0, -2.5 } // �p�^�[��2
{ -5.0, -5.0, -5.0, -5.0, 0.0 } // �p�^�[��3
{ -5.0, -5.0, -5.0, -5.0, 2.5 } // �p�^�[��4
{ -5.0, -5.0, -5.0, -5.0, 5.0 } // �p�^�[��5
{ -5.0, -5.0, -5.0, -2.5, -5.0 } // �p�^�[��6
{ -5.0, -5.0, -5.0, -2.5, -2.5 } // �p�^�[��7
{ -5.0, -5.0, -5.0, -2.5, 0.0 } // �p�^�[��8
						�E
						�E
						�E
{ 5.0, 5.0, 5.0, 5.0, 2.5 } // �p�^�[��3124
{ 5.0, 5.0, 5.0, 5.0, 5.0 } // �p�^�[��3125
*/
// ���񋓂���Ă���󋵂ɂȂ��Ă���͂�

// ���������p�^�[���̗񋓃��X�g�쐬
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

	make_que_wait(wait_que); // ���������̗񋓃��X�g����
	size_t counter = 0;

	std::mutex mtx;

	CsvWriter CW;
	CW.StandUp("Tester_");
	CW.fp7 << "�ڔ����Ԏ����ω���[sec],���c�n��[sec],�V��v��[sec],��X��[sec],���h[sec],�ϓd���o�͑��d�͗�[kWh]\n";

	std::vector<std::thread> TaskPool(0);
	TaskPool.reserve(NumThread);
	std::vector<SimulateEngine> se(NumThread);

	std::cout << "Starting Simulation \n";

	while (!wait_que.empty()) // ���X�g�ɃV�~�����[�V�������ׂ��������炵�p�^�[�����c���Ă���΁A���[�v���s
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
			TaskPool.at(i).join(); // ���������X���b�h�̏������ׂĂ������ł�������҂�
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


// main�֐����ǂ�
int main(void) {
	
	simulate();

	std::cout << "\n Finish Simulating\n Type to End\n";

	std::string waiter;
	std::cin >> waiter;

	return 0;

}
