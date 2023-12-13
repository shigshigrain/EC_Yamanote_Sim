#pragma once

#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <thread>
#include <future>
#include <memory>

#include "mySimulate.h"
#include "CsvWriter.hpp"


class SimulateEngine
{

private:
	/***  各行列宣言部  ***/
	/***  各行列宣言部  ***/
	std::unique_ptr<double[]> H; // = std::make_unique<double[]>(N_node * N_branch);		//Connection Matrix
	std::unique_ptr<double[]> H_tra; // = std::make_unique<double[]>(N_branch * N_node);	//Transpose matrix of H
	std::unique_ptr<double[]> Y; // = std::make_unique<double[]>(N_branch * N_branch);		//Conductance Matrix
	std::unique_ptr<double[]> A; // = std::make_unique<double[]>(N_node * N_node);			//A matrix(H*Y*H_tra)
	std::unique_ptr<double[]> TEMP; // = std::make_unique<double[]>(N_node * N_branch);	//Test
	std::unique_ptr<double[]> In; // = std::make_unique<double[]>(N_node);					//Node Current Matrix
	std::unique_ptr<double[]> In_cpy; // = std::make_unique<double[]>(N_node);				//Copy of Node Current Matrix
	std::unique_ptr<double[]> Vn; // = std::make_unique<double[]>(N_node);					//Node Voltage Matrix
	std::unique_ptr<double[]> Ib; // = std::make_unique<double[]>(N_branch);				//Branch Current Matrix
	std::unique_ptr<double[]> Vb; // = std::make_unique<double[]>(N_branch);
	std::unique_ptr<int[]> ipiv_Vn; // = std::make_unique<int[]>(N_node);
	//int i;
	//int j;

	int ERROR; // = 0;						//回路計算の再計算判定用フラグ（flagを全て加算してNUM_sub以上なら回路計算ループを抜ける。NUM_sub以下なら再計算。）
	int ERROR_REG; // = 0;
	int count_loop; // = 0;
	int count_loop_rec; // = 0;
	//	int test_flag = 0;

		/***  マトリックス計算関係  ***/
	const char trans = 'N'; //Normalな場合。これは規定値'N','T','C'から選択
	int K_node; // = N_node;
	int K_branch; // = N_branch;
	int ld_node; // = N_node;
	int ld_branch; // = N_branch;
	double alpha; // = 1.0;
	double beta; // = 0.0;
	//int* ipiv_Vn;

	int info_Vn; //;
	int nrhs; // = 1;
	int incx; // = 1;
	int incy; // = 1;
	/***  構造体宣言部  ***/
	std::unique_ptr<TRAIN[]> tra; // = std::make_unique<TRAIN[]>(NUM_tra);	//ポインタの配列として定義
	std::unique_ptr<SUB[]> sub; // = std::make_unique<SUB[]>(NUM_sub);
	std::unique_ptr<NODE[]> node; // = std::make_unique<NODE[]>(N_node);
	std::unique_ptr<NODE[]> node_order; // = std::make_unique<NODE[]>(N_node);
	std::unique_ptr<BRANCH[]> branch; // = std::make_unique<BRANCH[]>(N_branch);
	//std::unique_ptr<ESD[]> esd; // = std::make_unique<ESD[]>(NUM_esd);

	double timer;
	double minute;

	double sub_total;

public:
	SimulateEngine();
	bool Init();
private:
	void Make_train(TRAIN& tra, const INI_TRA& ini, int i);
	void Make_substation(SUB& sub, const INI_SUB& ini, int i);
	void Make_NODE_TRAIN(NODE& node, TRAIN& tra);
	void Make_NODE_SS(NODE& node, SUB& sub);
	void Error_Detection_SS(SUB& sub);
	void Error_Detection_REG(TRAIN& tra);
	void Calculate_R_run(TRAIN& tra);
	void Calculate_R_total(TRAIN& tra);
	void Traction_force(TRAIN& tra);
	void Regenerative_force(TRAIN& tra);
	double funcdelay(double aTs, double aTf, double ayold, double au);
	void d_current_power(TRAIN& tra);
	void d_current_regenerative(TRAIN& tra);
	void q_current_power(TRAIN& tra);
	void q_current_regenerative(TRAIN& tra);
	void d_q_voltage_calculation(TRAIN& tra);
	void Solve_next_state(TRAIN& tra);
	void Calculate_BEmax(TRAIN& tra);
	void change_direction(TRAIN& tra);
	void decide_final_station(TRAIN& tra, size_t i);
	void Run_pattern(TRAIN& tra, const std::vector<STATION>& sta, double t, std::vector<double>& WaitTimeIN, std::vector<double>& WaitTimeOUT);
	void Calculation_traction_force_brake_force(TRAIN& tra);
	void Calculation_traction_circuit(TRAIN& tra);
	void Calculation_power_motor(TRAIN& tra);
	void Calculation_power_train(TRAIN& tra);
	void Calculation_energy_train(TRAIN& tra);
	void Calculation_power_sub(SUB& sub);
	void Calculation_energy_sub(SUB& sub);
	void Sort_s(const std::unique_ptr<NODE[]>& data, int n);
	void Sort_m(const std::unique_ptr<NODE[]>& data, int n);
	void Initialize_BRANCH(BRANCH& branch);
	int Count_Trains_Direction(const std::unique_ptr<TRAIN[]>& tra, int direction);
	void Make_branch(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<NODE[]>& temp, const std::unique_ptr<TRAIN[]>& tra);
	void Initialize_matrix1(const std::unique_ptr<double[]>& Mtr, size_t M, size_t N);
	void Make_H_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Hm);
	void Make_Y_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Y);
	void Make_In_vector(const std::unique_ptr<NODE[]>& data, const std::unique_ptr<double[]>& In);
	void Transpose(const std::unique_ptr<double[]>& trans, const std::unique_ptr<double[]>& X, size_t row, size_t column);
public:
	int mySimulate(std::vector<double> WaitTimeIN, std::vector<double> WaitTimeOUT, CsvWriter& _cw, std::mutex& _mtx);

	~SimulateEngine();

};

// 非同期処理用非メンバ関数

bool ExeSimulate(SimulateEngine& SE, std::vector<double> _wtIN, std::vector<double> _wtOUT, CsvWriter& _cw, std::mutex& _mtx);