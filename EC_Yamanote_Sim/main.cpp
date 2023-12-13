
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

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w1�ł̒������Ԉꗗ^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w2^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w3^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w4^
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w5^
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w6�ȍ~���L�q�\�@�c�����̐��͉w�̐��ƈ�v^
};

const std::vector<std::vector<double>> ShiftWaitIN = {

	{-5.0, -4.0, 0.0}, // �ڔ�^
	{-5.0, -4.0, 0.0}, // ���c�n��^
	{-5.0, -4.0, 0.0}, // �V��v��^
	{-5.0, -4.0, 0.0}, // ��X��^
	{-5.0, -4.0, 0.0}, // ���h^
};

const std::vector<std::vector<double>> ShiftWaitOUT = {

	{-5.0, -4.0, 0.0}, // ���h^
	{-5.0, -4.0, 0.0}, // ��X��^
	{-5.0, -4.0, 0.0}, // �V��v��^
	{-5.0, -4.0, 0.0}, // ���c�n��^
	{-5.0, -4.0, 0.0}, // �ڔ�^
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
// ���񋓂���Ă���󋵂ɂȂ��Ă���͂�^

// ���������p�^�[���̗񋓃��X�g�쐬^
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

	make_que_wait(wait_IN, ShiftWaitIN); // ���������̗񋓃��X�g����^
	make_que_wait(wait_OUT, ShiftWaitOUT); // ���������̗񋓃��X�g����^
	size_t counter = 0;

	std::mutex mtx;

	CsvWriter CW;
	CW.StandUp("Tester_");
	CW.fp7 << "��:�ڔ���Ԏ����ω�[s],��:���c�n��[s],��:�V��v��[s],��:��X��[s],��:���h[s],"
		   << "�O:���h��Ԏ����ω�[s],�O:��X��[s],�O:�V��v��[s],�O:���c�n��[s],�O:�ڔ�[s],"
		   << "�ϓd���o�͑��d�͗�[kWh]\n";

	std::vector<std::thread> TaskPool(0);
	TaskPool.reserve(NumThread);
	std::vector<SimulateEngine> se(NumThread);

	std::cout << "Starting Simulation \n";

	// �����E�O��肻�ꂼ��̃��[�v�Ŏ���^
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
				tp.join(); // ���������X���b�h�̏������ׂĂ������ł�������҂�
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


	//while (!wait_que.empty()) // ���X�g�ɃV�~�����[�V�������ׂ��������炵�p�^�[�����c���Ă���΁A���[�v���s
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
	//		TaskPool.at(i).join(); // ���������X���b�h�̏������ׂĂ������ł�������҂�
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
