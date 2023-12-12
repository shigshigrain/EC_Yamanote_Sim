#pragma once

/********************************************************************
Simulation of train running and feeder circuit operation

ver13_2
1 train, 4 substation

ver14
Simulate Takasaki line
(With Diagram function)

ver15
Simulate Takasaki line
(With Diagram function)
With Wayside ESD
���ёւ��̃A���S���Y�����������iOFF��SS��Ђ��[����E�[�Ɉړ�����������j

ver16
Simulate Takasaki line
(With Diagram function)
With Wayside ESD & regenerative INV

ver16_1
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV

ver16_2
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV

ver16_3
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV
With Delay of departure

ver16_4
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD(time limitation included) & regenerative INV
With Delay of departure
With INV-ESD interruction investigated


ver17
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD(time limitation included) & regenerative INV
With Delay of departure
With INV-ESD interruction investigated
�p�����[�^�����ύX�@�\�ǉ��H�i�j
ESD: V_charge_ESD, V_discharge_ESD
INV: V_charge_INV (V_fin)

ver18�@���ѐ搶�̂��̂���������
���ѐ搶�̂͌��X�C�񐶃u���[�L��b_reg_c[kN]�͌v�Z���ĂȂ����ǁC�񐶃u���[�L��b_reg[N]�͌v�Z���ĂĂ�����񐶃G�l���M�[�̌v�Z�ɗp���Ă���

�m�b�`5�͍s�� v1_P=40,v2_P=63.82549701  i1d=85.7418508,iq=165.6086836,ws=11.28173426,Fmax=13.90011356

�񐶎�v1_R=50.52,i1d=103.2188528,i1q=-137.5678441,ws=-7.784729348
Bmax=13.90011356


�m�b�`4�񐶎�v1_R=56.31,i1d=102.7443702,i1q=-111.2973885,ws=-6.32721445

�m�b�`3
�񐶎�i1d=102.2117125,i1q=-84.83142392,ws=-4.84776665

d���d����q���d���v�Z��ǉ�����
�������烂�[�^�p���[�����߂邱�ƂŁC�C���o�[�^�����𖳎����C�C���o�[�^�d�����v�Z
���d�d���ϓ����C�͍s���̈����͂��ς��悤�ɂ���
�񐶎��̉񐶎��̃u���[�L�͂�ύX(�y���׉񐶐���ɂ�萫�\���ω�)

ver19�@2022.03.27
���c�}����Ƃ̋��������ŏ��c�}�������ł̑��s�V�~�����[�V���������邱�ƂɂȂ���(����)
�Ƃ肠�����V�S�����u�`���ؓc�Ԃ̂��̂ɓK�p�����Ă݂�

********************************************************************/


//char dir[]= "C:\Users\�G��\Desktop\��������\JRE\����\2018_10\\";


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <math.h>
#include "mkl.h"
#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include "CsvWriter.hpp"


using namespace std;

/***  �V�~�����[�V�������ݎ���  ***/
//#define dt 0.01				//���ݎ���:10ms

#define dt 0.001			//���ݎ���:1ms

//#define dt 0.0005			//���ݎ���:500��s

//#define dt 0.005
//#define dt 0.0001			//���ݎ���:100��s	//�����蒷��ƃG���[���N����(CFC������������)
//#define dt 0.00005		//���ݎ���:50��s
//#define dt 0.00001		//���ݎ���:10��s

//double t, minute;

//#define NUM_tra_UP 1								//Number of UP trains(+1)
//#define NUM_tre_DOWN 1								//Number of DOWN trains(-1)
//#define NUM_tra (NUM_tra_UP+NUM_tre_DOWN)			//Number of trains

#define NUM_tra 2										//Number of trains
#define NUM_sub 3										//Number of substations
#define NUM_esd 0										//Number of ESDs
#define NUM_inv 0										//Number of regenerative inverters
#define NUM (NUM_tra + NUM_sub + NUM_esd + NUM_inv - 1)						//Number of elements
#define MAX_tra 24
#define N_node (NUM_tra+NUM_sub)					//Number of nodes
#define N_branch (NUM_sub + 2*(NUM_sub-1)+NUM_tra)	//Number of branches

#define NUM_station 8
#define NUM_final_station 8			//�e�w��ԁ@(�w�� - 1 �������H)

#define NUM_station2 5				//�}�s(�w�� - 1 �������H)
#define NUM_final_station2 5

#define Ratio_E	3600.0*1000.0							//kwh��J�ւ̕ϊ���

#define T_charge 50.0
#define T_discharge 40.0
#define T_ope 180.0



/** �o�̓f�[�^�t�@�C�� **/
				//��Ԏ����ύX���E�ϓd���o�͑��d�͗��o�͗p�t�@�C��

#define Consider_loop_flag 0 //1�Ȃ烋�[�v���Ȃ�

//#define Vss 1650.0						//�����׎��ϓd���o�͓d��[V]
//#define Vss 1600.0						//�����׎��ϓd���o�͓d��[V]
#define Vss1 1506.0						//�����׎��ϓd���o�͓d��[V]
#define Vss2 1532.9						//�����׎��ϓd���o�͓d��[V]
#define Vss3 1497.6						//�����׎��ϓd���o�͓d��[V]
//#define Vss 1500.0						//�����׎��ϓd���o�͓d��[V]
//#define Iss_rated 8000.0						//�ϓd����i�d��[A]
#define Iss_rated1 8000.0						//�ϓd����i�d��[A]
#define Iss_rated2 4000.0						//�ϓd����i�d��[A]
#define Iss_rated3 4000.0						//�ϓd����i�d��[A]

#define ess1 5.78						//�ϓd���o�͓d���ϓ���[%]
//#define ess1 5.89						//�ϓd���o�͓d���ϓ���[%]
#define ess2 5.59						//�ϓd���o�͓d���ϓ���[%]
#define ess3 7.89						//�ϓd���o�͓d���ϓ���[%]

#define R	0.0360959					//���d��R[��/km]
#define R_StoN_UP	0.03407					//���d��R[��/km]
#define R_StoN_DOWN	0.03423					//���d��R[��/km]
#define R_NtoK_UP	0.04018					//���d��R[��/km]
#define R_NtoK_DOWN	0.04061					//���d��R[��/km]

#define C_FC (10.05 / 1000.0)		//�ԗ�FC�e��[F](5M��)
#define L_FL 0.0080				//FL����[H](5M��)
#define R_FL 0.0348				//FL������R[��](5M��)
//#define C_FC 0.015		//�ԗ�FC�e��[F](5M��)
//#define L_FL 0.019				//FL����[H](5M��)
//#define R_FL 0.100				//FL������R[��](5M��)


#define Vclim	1720.0
#define Vcmax	(Vclim + 30.0)

#define g 9.80665						//�d�͉����x[m/s/s]

#define speed_limit 90.0				//�������x[km/h]
#define speed_min   50.0                //�ĉ������x[km/h]

#define stoptime 40.0					//��Ԏ���[s]

//#define P_SIV 142.0*3.0*1000.0/10.0				//��@�d��[W]�i1��������j
//#define P_SIV 300.0*1000.0/10.0				//��@�d��[W]�i1��������j
#define P_SIV 20.0*1000.0/10.0
//#define P_SIV 0.0

#define P_KARAKIDA 24.0*445.0*2.0*5.0	//�o�͓d���~�o�͓d���~�����~�Ґ���
#define P_SHINYURIGAOKA 1000.0*1500.0

#define count_lap 1

#define rownum (int)(speed_limit+15)*20


/******************* �o�b�e���[�Ɋւ���ݒ� ********************/
/*LIM50EN-8*/
#define NoS				4.0							//����
#define NoP				36.0						//����
#define Vbat_mod_rated	173.0						//���W���[����i�d��[V]
#define Ibat_mod_rated	200.0						//���W���[����i�d��[A]
#define	Ebat_mod_c		0.951						//1���W���[��������G�l���M�[�e��[kWh]
#define Rbat_mod		0.5							//1���W���[�������������R[��]
#define C_rate_mod		2.0							//1���W���[��������C���[�g(�[�d�F2.5C�A���d�F6C)

#define Vbat_rated		Vbat_mod_rated*NoS			//�~�d�f�q�d��[V]
#define Eesd_cap_c		Ebat_mod_c*NoS*NoP			//�~�d�f�q�G�l���M�[�e��[kWh]
#define Eesd_cap		Eesd_cap_c*Ratio_E			//�~�d�f�q�G�l���M�[�e��[J]
#define Rbat			Rbat_mod*NoS/NoP			//�~�d�f�q������R[��]
#define C_rate			C_rate_mod*NoP				//�~�d�f�qC���[�g
/*************************************************************/



/****** ����萔�֌W�i�d������j ******/
#define Vc_min 1400.0
#define Vc_max Vcmax

#define lim_esd_c 1.0		//[kWh]

#define T_control 0.00001		//���䉉�Z����[s]�i���ݎ��ԂƓ����j
/**************************************/


/****** ����~�d���u����p�����[�^ ******/
#define V_charge_ESD 1715.0			//�[�d�J�n�d��
#define V_discharge_ESD 1665.0		//���d�J�n�d��
/**************************************************/


/****** ����񐶃C���o�[�^����p�����[�^ ******/
#define V_charge_INV 1675.0			//����J�n�d��
#define V_fin 1673.0				//����I���d��
/**************************************************/


/****** �e����� ******/
#define e_gear 0.98		//�M�A����
#define e_motor 0.92	//���[�^�[����
#define e_inv 0.95		//�C���o�[�^����
/**********************/


/////�X�C�b�`���O�f�q�̑����v�Z�p�p�����[�^/////////////////
/*** CM1200HC-66H ***/
#define t_on 1.6/1000.0/1000.0
#define t_off 2.5/1000.0/1000.0
#define Vce 3.6
#define Vf 2.7
#define Irp 1200.0
#define t_rr 1.4/1000.0/1000.0
#define fc 1000.0


//double t_out;		//�f�[�^�Ԉ����֌W
//double t_out_old;
//double t_out2;		//�V�~�����[�V�����i���\���֌W
//double t_out2_old;
//double t_ope;
//
//int loop_endflag = 0;
//int sim_endflag = 0;


typedef struct {
	/**  �萔  **/
	int name_train;			//��Ԃ̖��O
	double mass;			//��ԑ��d��[kg]
	double mass_c;			//��ԑ��d��[t]
	double mass_M;			//M�ԏd��[kg]
	double mass_T;			//T�ԏd��[kg]
	double mass_P;			//�׏d�d��[kg]
	double num;				//����
	double num_M;			//M�Ԃ̗���
	double num_T;			//T�Ԃ̗���
	double num_mot;			//���[�^��
	double v0_P;           //(�͍s��)V/f�I�[���x
	double v1_P;			//(�͍s��)��o�͗̈�ɓ��鑬�x(5�m�b�`)
	double v2_P;			//(�͍s��)�����̈�ɓ��鑬�x(5�m�b�`)
	double v0_R;           //(�񐶎�)V/f�I�[���x
	double v1_R;			//(�񐶎�)��o�͗̈�ɓ��鑬�x(5�m�b�`)
	double v2_R;			//(�񐶎�)�����̈�ɓ��鑬�x(5�m�b�`)
	double Fmax;			//��g���N�̈�̈�����(�͍s��)
	double Bmax;			//��g���N�̈�̈�����(�񐶎�)
	int  direction;			//���E���蔻��(1:��{�ː�����ʁC-1:����ˑ�{����)
	double UP_DOWN;
	double B;				//��p�ő匸���x[m/s/s]
	double Bref;			//�����x�w��[m/s/s]
	double BEref;			//�u���[�L�͎w��[N]
	double BEmax;			//�ő�u���[�L��[N]
	double Cfc;//FC�e��
	double Lfl;//FL�C���_�N�^���X
	double Rfl;//Fl��R

	double theta;
	double theta_old;

	double notch;				//�m�b�`����(Powering: 1�`4)
	double brake;				//�m�b�`����(Regenerating: -1�`-7)
	int Type;				//��Ԏ�ʁC0:�e�w��ԁC1:�}�s
	int t_count;			//���^�]�ԉ^�]���@
	int c_count;			//�葬�^�]�J�E���g
	double con_speed;		//�葬�^�]���x

	/**  �ϐ�  **/
	int accelflag;			//accelflag=1:�͍s, accelflag=2:�čs, accelflag=3:����, accelflag=4:���, accelflag=5:���x����
	int brakeflag;			//brakeflag=1:�w�Ɏ~�܂邽�߂Ɍ�����, brakeflag=0:����ȊO
	int laststopflag;		//laststopflag=1:�I�_�ɓ���, laststopflag=0:����ȊO
	double X_recentstop;		//���O�̒�Ԉʒu[m]
	double X_nextstop;		//���̒�Ԉʒu[m]
	double X_brake;			//�u���[�L�J�n�ʒu[m]
	double x;				//�ʒu[m]
	double x_c;				//�ʒu[km]
	double v;				//���x[m/s]
	double v_c;				//���x[km/h]
	double v_new;
	double a;				//�����x[m/s/s]
	double a_c;				//�����x[km/h/s]
	double Fmot;			//������[N]
	double Fmot_c;			//������[kN]
	double Fmot_regM;		//�ő�񐶃u���[�L��[N]
	double Fmot_regM_c;		//�ő�񐶃u���[�L��[kN]
	double Ftot;			//�g�[�^��������[N]
	double Ftot_c;			//�g�[�^��������[kN]
	double b_reg;			//�񐶃u���[�L��[N]
	double b_reg_c;			//�񐶃u���[�L��[kN]
	double b_reg_judge;		//�y���׉񐶐���p�^�[���𖞂����Ă��邩�̔���p
	double judgement;		//�J��Ԃ����ۂ�
	double b_air;			//�@�B�u���[�L��[N]
	double b_air_c;			//�@�B�u���[�L��[kN]
	double b_reg_loss;			//�񐶃u���[�L��[N]
	double b_reg_loss_c;			//�񐶃u���[�L��[kN]

	double i1d_P; //�͍s���C��g���N�̈��d���d��
	double i1q_P; //�͍s���C��g���N�̈��q���d��
	double i1d_R; //�񐶎��C��g���N�̈��d���d��
	double i1q_R; //�񐶎��C��g���N�̈��q���d��
	double i1q_R_1;//�񐶎��C��o�͗̈悪�Ȃ������ōő�ƂȂ�d��
	double R1;//1����R[��]
	double R2;//2����R[��]
	double L1;//1�����ȃC���_�N�^���X[H]
	double L2;//2�����ȃC���_�N�^���X[H]
	double MORE1;//�ꎟ�R��C���_�N�^���X[H]	
	double MORE2;//�񎟘R��C���_�N�^���X[H]
	double g0;//�㎥�R���_�N�^���X[S]
	double M_IM;//���݃C���_�N�^���X[H]
	double SIGMA;//�R��W��
	double Gr;//�M�A��
	double rd;//�ԗ֔��a[m]
	double POLE;//�ɑΐ�
	double Ts;//��������
	double Tf;//�����萔

	double id_mot_P;       //�͍s����d���d��[A]
	double iq_mot_P;       //�񐶎���q���d��[A]
	double id_mot_R;       //�͍s����d���d��[A]
	double iq_mot_R;       //�񐶎���q���d��[A]
	double idtot;          //�͍s�E�čs�E�񐶎���d���d��[A]
	double iqtot;          //�͍s�E�čs�E�񐶎���q���d��[A]
	double idtotF;          //�͍s�E�čs�E�񐶎���d���d��[A]
	double iqtotF;          //�͍s�E�čs�E�񐶎���q���d��[A]
	double iq_reg;         //�񐶎���q���d��(�y���׉񐶐��䍞��)[A]
	double iq_regF;         //�񐶎���q���d��(�y���׉񐶐��䍞��)[A]

	double wr;//�p���g��[rad/s]
	double ws;//���ׂ�p���g��[rad/s]
	double we;//�C���o�[�^�p���g��[rad/s]
	double v1d;//d���d��[V]
	double v1q;//q���d��[V]
	double Em;//�C���o�[�^�o�͓d��[V]

	double Pmot;			//���[�^�p���[[W]
	double Pmot_air;		//�@�B�u���[�L�p���[[W]
	double Pmot_c;			//���[�^�p���[[kW]
	double Pmot_air_c;		//�@�B�u���[�L�p���[[kW]

	double Pmot_reg_loss;
	double Pmot_reg_loss_c;

	double Ptot;			//�ԗ��p���[[W](���[�^�p���[+�@�B�u���[�L��)
	double Ptot_c;			//�ԗ��p���[[kW](���[�^�p���[+�@�B�u���[�L��)

	double Pres;			//���s��R�ł̑���[W]
	double Pres_c;			//���s��R�ł̑���[kW]

	double Pinv_on;			//�X�C�b�`���O�f�q�̃I������[W]
	double Pinv_diode;		//�_�C�I�[�h�̓��ʑ���[W]
	double Pinv_sw;			//�X�C�b�`���O����[W]
	double Pinv_rec;		//�_�C�I�[�h�̃��J�o���[����[W]

	double Ploss_inv;		//�C���o�[�^����[W]
	double Ploss_mot;		//���[�^����[W]
	double Ploss_fe;		//���[�^�S��[W]
	double Ploss_fl;		//FL����[W]
	double Ploss_all;		//�����H����[W]

	double Ploss_inv_c;
	double Ploss_mot_c;
	double Ploss_fl_c;
	double Ploss_all_c;

	double Pvh;				//�ԗ����o�̓p���[[W] (Pvh_in + Pvh_out)
	double Pvh_in;			//�ԗ����̓p���[[W]
	double Pvh_out;			//�ԗ��o�̓p���[[W]
	double Pvh_st;			//�ԗ����̓p���[(�čs�A��Ԏ�)[W]
	double Pvh_c;			//�ԗ����o�̓p���[[kW]
	double Pvh_in_c;		//�ԗ����̓p���[[kW]
	double Pvh_out_c;		//�ԗ��o�̓p���[[kW]
	double Pvh_st_c;			//�ԗ����̓p���[(�čs�A��Ԏ�)[kW]

	double vfc;				//FC�d��[V]
	double iinv;				//�C���o�[�^�d��[A]
	double ifc;				//FC�d��[A]
	double ifl;				//FL�d��[A]
	double vfc_old;				//FC�d�������O[A]
	double isiv;			//FL�d��[A]
	double vp;				//�ԗ��p���^�_�d��[V]
	double vp_old;				//�ԗ��p���^�_�d��[V]
	double R_run;			//���s��R[N/t]
	double R_grade;			//���z��R[N/t]
	double R_curve;			//�Ȑ���R[N/t]
	double R_total;			//�S���s��R[N]

	/** �G�l���M�[�v�Z�֌W **/
	double Emot_pow;		//���[�^�͍s�G�l���M�[[J]
	double Emot_reg;		//���[�^�񐶃G�l���M�[[J]
	double Emot_air;		//�@�B�u���[�L�G�l���M�[[J]
	double Emot_reg_loss;		//���[�^�񐶍i�荞�݃G�l���M�[[J]

	double Emot_pow_c;		//���[�^�͍s�G�l���M�[[kWh]
	double Emot_reg_c;		//���[�^�񐶃G�l���M�[[kWh]
	double Emot_air_c;		//�@�B�u���[�L�G�l���M�[[kWh]
	double Emot_reg_loss_c;		//���[�^�񐶍i�荞�݃G�l���M�[[kWh]

	double Eloss_inv_pow;	//�C���o�[�^�����G�l���M�[(�͍s��)[J]
	double Eloss_mot_pow;	//���[�^�����G�l���M�[(�͍s��)[J]
	double Eloss_fl_pow;	//FL�����G�l���M�[(�͍s��)[J]
	double Eloss_all_pow;	//�����H�����G�l���M�[(�͍s��)[J]

	double Eloss_inv_reg;	//�C���o�[�^����(�񐶎�)[J]
	double Eloss_mot_reg;	//���[�^����(�񐶎�)[J]
	double Eloss_fl_reg;	//FL����(�񐶎�)[J]
	double Eloss_all_reg;	//���H����(�񐶎�)[J]

	double Eres_pow;		//���s��R�ł̑���(�͍s��)[J]
	double Eres_reg;		//���s��R�ł̑���(�񐶎�)[J]
	double Eres_coa;		//���s��R�ł̑���(�čs��)[J]

	double Esiv_pow;	//SIV�̏���G�l���M�[(�͍s��)[J]
	double Esiv_reg;	//SIV�̏���G�l���M�[(�񐶎�)[J]
	double Esiv_coa;	//SIV�̏���G�l���M�[(�čs��)[J]
	double Esiv_stp;	//SIV�̏���G�l���M�[(��Ԏ�)[J]

	double Eloss_inv_pow_c;
	double Eloss_mot_pow_c;
	double Eloss_fl_pow_c;
	double Eloss_all_pow_c;

	double Eloss_inv_reg_c;
	double Eloss_mot_reg_c;
	double Eloss_fl_reg_c;
	double Eloss_all_reg_c;

	double Eres_pow_c;
	double Eres_reg_c;
	double Eres_coa_c;

	double Esiv_pow_c;		//SIV�̏���G�l���M�[(�͍s��)[kWh]
	double Esiv_reg_c;		//SIV�̏���G�l���M�[(�񐶎�)[kWh]
	double Esiv_coa_c;		//SIV�̏���G�l���M�[(�čs��)[kWh]
	double Esiv_stp_c;		//SIV�̏���G�l���M�[(��Ԏ�)[kWh]

	double Evh;				//��ԓ��o�̓G�l���M�[[J]
	double Evh_in;			//��ԓ��̓G�l���M�[[J]
	double Evh_out;			//��ԏo�̓G�l���M�[[J]
	double Evh_st;          //��ԓ��̓G�l���M�[(�čs�A��Ԏ�)[J]

	double Evh_c;			//��ԓ��o�̓G�l���M�[[kWh]
	double Evh_in_c;		//��ԓ��̓G�l���M�[[kWh]
	double Evh_out_c;		//��ԏo�̓G�l���M�[[kWh]
	double Evh_st_c;		//��ԓ��̓G�l���M�[(�čs�A��Ԏ�)[kWh]

	/** ��Ԓ�~�ʒu�v�Z�֌W **/
	int flag_station;			//��ԉw��
	int startflag;				//���Ԏw��
	/**  ��Ԓ�Ԏ��Ԋ֌W  **/
	double t_stop;				//
	double t_stop_old;			//�w���B����[s]
	double t_wait;
	double route[rownum][2];	//�u���[�L���[�g{�ʒu,���x}
	int diaflag;				//�^�s����
	int lap;					//����
	int lapmax;					//�ڕW����
	int reg_flag;				//�y���׉񐶐���t���O

	double X_nextstop_old;		//���̒�Ԉʒu[m]
	double Speedlimit;			//�������x[km/h]
	double Reaccelspeed;			//�������x[km/h]

	double T_delay;				//���ԃ^�C�~���O�x��l��[s]
	double Vst;					//���x�����͓����x�[�X�d��

} TRAIN;

typedef struct {
	/** �萔 **/
	int name_SS;		//�ϓd���̖��O
	double Vss_0;		//�����׎��ϓd���o�͓d��
	double Iss_0;		//��i�o�͓d��
	double e_ss;		//�ϓd���o�͓d���ϓ���
	double Xss;			//�ϓd���ʒu
	double Jss_ini;		//�ϓd�������d����[�����l]
	double Rss_ini;			//�ϓd�������o�͒�R[�����l]

	/** �ϐ� **/
	int diode;			//�ϓd���_�C�I�[�hON/OFF����t���O(1:ON 0:OFF)
	double Rss;			//�ϓd�������o�͒�R[��]
	double Jss;			//�ϓd�������d����[A](�X�V�p)

	double vss;			//�ϓd���o�͓d��[V]
	double vss_e;		//�ϓd���o�͓d��[V]
	double iss;			//�ϓd���o�͓d��[A]
	double Pss;			//�ϓd���o�̓p���[[W]
	double Pss_c;		//�ϓd���o�̓p���[[kW]
	double vout;		//�ϓd������d��[V]
	int flag;			//�t���O

	/** �G�l���M�[�v�Z�֌W **/
	double Ess;			//�ϓd���o�̓G�l���M�[[J]
	double Ess_c;		//�ϓd���o�̓G�l���M�[[kWh]
	double Wss;			//�ϓd������[J]
	double Wss_c;		//�ϓd������[kWh]
} SUB;

typedef struct {
	/**  �萔  **/
	char name_station[20];	//�w�̖��O
	double Xs;	//����
} STATION;

typedef struct {
	/** �萔 **/
	int name_ESD;		//ESD�̖��O
	double Cfc;			//FC�Ód�e��[F]
	double RL;			//���A�N�g��������R[��]
	double RFL;			//FL������R[��]
	double SOC_min;		//SOC�����l[%]
	double SOC_max;		//SOC����l[%]

	/** ����萔 **/
	/*�d������֌W*/
	double kp_v_ch;			//�d������n���Q�C��(�[�d��)
	double ki_v_ch;			//�d������n�ϕ��Q�C��(�[�d��)
	double kp_v_disch;		//�d������n���Q�C��(���d��)
	double ki_v_disch;		//�d������n�ϕ��Q�C��(���d��)
	double kp_v;			//�d������n���Q�C��
	double ki_v;			//�d������n�ϕ��Q�C��
	double P_v;				//�d���������o��
	double I_v;				//�d������ϕ���o��

	double Eesd_min;				//ESD�G�l���M�[�����l[J]
	double Eesd_max;				//ESD�G�l���M�[����l[J]
	double Eesd_min_c;				//ESD�G�l���M�[�����l[kWh]
	double Eesd_max_c;				//ESD�G�l���M�[����l[kWh]
	double Eesd_minlim0;			//[kJ]�@�G�l���M�[���������������l
	double Eesd_maxlim0;			//[kJ]�@�G�l���M�[��������������l
	double Eesd_minlim0_c;			//[kWh]�@�G�l���M�[���������������l
	double Eesd_maxlim0_c;			//[kWh]�@�G�l���M�[��������������l

	double Vc_mid;			//���ԓd��
	double Eesd_minlim;		//[J]�@�G�l���M�[���������������l
	double Eesd_maxlim;		//[J]�@�G�l���M�[��������������l
	double Eesd_minlim_c;	//[kWh]�@�G�l���M�[���������������l
	double Eesd_maxlim_c;	//[kWh]�@�G�l���M�[��������������l


	/*�d������֌W*/
	double kp_i;			//�d������n���Q�C��
	double ki_i;			//�d������n�ϕ��Q�C��
	double P_i;				//�d������n����o��
	double I_i;				//�d������n�ϕ���o��

	double V_charge;		//[V]�@���d�J�n�d��
	double V_discharge;		//[V]�@�[�d�J�n�d��
	double SOC_ref;			//SOC�w�ߒl�i�~�d�ʊǗ��p�j
	double Iesd_mid;		//SOC�Ǘ��̂��߂̒����[���d�p



	double kel;
	double keh;
	double iesd_ref;		//�������d���w�ߒl
	double v0_ref;			//�������d���w�ߒl
	double vp_ref;			//�ː��d���w�ߒl
	double a;				//�ʗ���(vesd/v0)

	double Pesd_charge_lim;			//[W]�@�[�d�p���[�ő�l
	double Pesd_charge_lim_c;		//[kW]�@�[�d�p���[�ő�l
	double Pesd_discharge_lim;		//[W]�@���d�p���[�ő�l
	double Pesd_discharge_lim_c;	//[kW]�@���d�p���[�ő�l
	double Iesd_max;		//��iESD�d��[A]
	double Imax;			//�ő�[���d�d��[A]


	/** �ϐ� **/
	double x;			//ESD�ݒu�ʒu[m]
	double vp;			//�p���^�_�d��[V]
	double idc;			//�������d��[A]
	double vc;			//FC�d��[V]
	double ic;			//FC�d��[A]
	double ich;			//�`���b�p�o�͓d��[A]
	double v0;			//IGBT�R���N�^-�G�~�b�^�ԓd��
	double vesd;		//ESD�d��[V]
	double iesd;		//ESD�d��[A]
	double Pesd;		//ESD�p���[[W]
	double Pesd_c;		//ESD�p���[[kW]
	double Pesd_loss;	//ESD�p���[[W]
	double Pesd_loss_c;	//ESD�p���[[kW]
	double Eesd;		//ESD�G�l���M�[[J]
	double Eesd_c;		//ESD�G�l���M�[[kWh]
	double Eesd_loss;	//ESD�����G�l���M�[[J]
	double Eesd_loss_c;	//ESD�����G�l���M�[[kWh]

	double Eesd_charge;			//ESD�[�d�G�l���M�[[J]
	double Eesd_charge_c;		//ESD�[�d�G�l���M�[[kWh]
	double Eesd_discharge;		//ESD���d�G�l���M�[[J]
	double Eesd_discharge_c;	//ESD���d�G�l���M�[[kWh]

	double SOC;			//[%]

	double t_charge;			//�[���d���Ԋ֌W
	double t_discharge;
	double t_ope;
} ESD;

typedef struct {
	/** �萔 **/
	int name_INV;		//��INV�̖��O
	double Cfc;			//FC�Ód�e��[F]
	double RL;			//���A�N�g��������R[��]
	double RFL;			//FL������R[��]

	/** ����萔 **/
	/*�d������֌W*/
	double kp_v;			//�d������n���Q�C��
	double ki_v;			//�d������n�ϕ��Q�C��
	double P_v;				//�d���������o��
	double I_v;				//�d������ϕ���o��

	/*�d������֌W*/
	double kp_i;			//�d������n���Q�C��
	double ki_i;			//�d������n�ϕ��Q�C��
	double P_i;				//�d������n����o��
	double I_i;				//�d������n�ϕ���o��

	double Vp_start;		//[V]�@����J�n�d��
	double Vp_fin;		//[V]�@����I���d��


	double iinv_ref;		//�������d���w�ߒl
	double vp_ref;			//�ː��d���w�ߒl

	double Pinv_lim;		//[W]�@�[���d�p���[�ő�l
	double Pinv_lim_c;		//[kW]�@�[���d�p���[�ő�l
	double Iinv_max;		//��INV �ő�d��
	double Iinv_rated;		//��INV ��i�d��


	/** �ϐ� **/
	double x;			//ESD�ݒu�ʒu[m]
	double vp;			//�p���^�_�d��[V]
	double idc;			//�������d��[A]
	double vc;			//FC�d��[V]
	double ic;			//FC�d��[A]
	double iinv;		//��INV �o�͓d��[A]
	double Pinv;		//��INV �񐶃p���[[W]
	double Pinv_c;		//��INV �񐶃p���[[kW]
	double Einv;		//��INV �񐶃G�l���M�[[J]
	double Einv_c;		//��INV �񐶃G�l���M�[[kWh]
} INV;

struct INI_TRA {
	double mass_M;
	double mass_T;
	double num_M;
	double num_T;
	double x;
	double v;
	int direction;
	double T_delay;
	int Type;
	int lapmax;
};

struct INI_SUB {
	double Xss;
	double e_ss;
	double Iss;
	double vss;
	double iss;
	int diode;
};


typedef struct {
	int Number;			//�m�[�h�ԍ�
	double X;			//�ʒu[m]�i��_����̋����j
	double I;			//�m�[�h�d��[A]
	double V;			//�m�[�h�d��[V]
	double r;			//�m�[�h��R�l[��](Only SS)
	int flag;		//�m�[�h�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
} NODE;

typedef struct {
	double X_start;		//branch�n�[[m]�i��_����̋����j
	double X_end;		//branch�I�[[m]�i��_����̋����j
	double r;			//Resistance [Ohm]
	int flag;		//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
	int Node_pos;	//���������ɐڑ�����Ă���m�[�h��
	int Node_neg;	//���������ɐڑ�����Ă���m�[�h��
} BRANCH;

/*** ��ԍ\���̔z�񏉊����p�̍\���̔z��錾�� ***/
	/*�I�t�s�[�N�_�C��*/
	//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���, ��Ԏ��
	//�R�����14�F45-15:30����ɂ���
	//�����F�r�܁��V�h���a�J	

const std::vector<INI_TRA> ini_tra = {
	//{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 0.0, 0.0, 1, 1080.0, 0, 0},				//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
	{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 0.0, 0.0, 1, 1380.0, 0, 0},				//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
	//{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 8200.0, 0.0, -1, 1380.0, 0, 0},				//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
	{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 8200.0, 0.0, -1, 1440.0, 0, 0},		//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
	//{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 0.0, 0.0, 1, 1440.0, 0, 0},		//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���

};

/*** �ϓd���\���̔z�񏉊����p�̍\���̔z��錾 ***/
const std::vector<INI_SUB> ini_sub = {
	{-300.0, ess1, Iss_rated1, Vss1, 0.0, 1},						//[�V�S�����u]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d��, �ː��d���C�_�C�I�[�h
	{5100.0, ess2, Iss_rated2, Vss2, 0.0, 1},						//[�i�R(��)]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d��, �ː��d���C�_�C�I�[�h
	{8500.0, ess3, Iss_rated3, Vss3, 0.0, 1},					//[���ؓc]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d���C�ː��d���C�_�C�I�[�h
};

/***  �w�̍\���̔z��錾�� (�V�S�����u�w�) ***/
const STATION sta[] = {	//����e�w
	{"basement1",	-100.0},			//�w��, �w�ʒu	0
	{"Ikebukuro",  	0.0},				//�w��, �w�ʒu	1
	{"Mejiro",	    1200.0},	    	//�w��, �w�ʒu	2
	{"Takadanobaba",  2100.0},			//�w��, �w�ʒu	3
	{"Shin-Okubo",	3500.0},			//�w��, �w�ʒu	4
	{"Shinjuku",	4800.0},			//�w��, �w�ʒu	5
	{"Yoyogi",	5500.0},		    	//�w��, �w�ʒu	6
	{"Harajuku",	7000.0},			//�w��, �w�ʒu	7
	{"Shibuya",	8200.0},			    //�w��, �w�ʒu	8
	{"basement2",	8300.0},			//�w��, �w�ʒu	9
};



//void Make_file7(FILE** fp)
//{
//	time_t start;
//	errno_t error;
//	struct tm today;
//	char filename[128];
//
//	time(&start);
//	error = localtime_s(&today, &start);
//
//	sprintf_s(filename, sizeof(filename), "Result\\2���%d��%d��%d��%d��%d�b.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//���O�̍쐬
//	error = fopen_s(fp, filename, "w");																														//�t�@�C���I�[�v��
//}
//
//
//void Make_train(TRAIN& tra, const INI_TRA* ini, int i)
//{
//	int m;
//	int n;
//
//	//*tra = (TRAIN*)malloc(sizeof(TRAIN));
//	//if (*tra == NULL)
//	//{
//	//	//�G���[����
//	//}
//
//	tra.name_train = i;
//	tra.Type = ini->Type;
//	tra.notch = 3.0;
//	tra.brake = 5.0;
//
//	tra.mass_M = ini->mass_M * ini->num_M;			//[kg]�����p(100%���)
//	tra.mass_T = ini->mass_T * ini->num_T;			//[kg]�����p(100%���)
//	tra.mass_P = (ini->num_M * 152.0 + (ini->num_T - 2.0) * 152.0 + 2.0 * 144.0) * 55.0 / 2.0;
//
//	if (tra.Type == 2) tra.mass_P = 0.0;
//
//	tra.mass = tra.mass_M + tra.mass_T + tra.mass_P + 15.8 * 1000.0 * (ini->num_M + ini->num_T) / 10.0;				//[kg](250%���)
//	tra.mass_c = tra.mass / 1000.0;										//[t](250%���)
//
//	tra.num = ini->num_M + ini->num_T;
//	tra.num_mot = ini->num_M * 4.0;
//
//	tra.v0_P = 38.0 / 3.6;        //[m/s]
//	tra.v1_P = 44.0 / 3.6;		//[m/s]
//	tra.v2_P = 46.0 / 3.6;		//[m/s]
//	tra.v0_R = 52.0 / 3.6;        //[m/s]
//	tra.v1_R = 86.0 / 3.6;		//[m/s]
//	tra.v2_R = 85.0 / 3.6;		//[m/s]
//
//	if (tra.Type == 2) {
//		tra.v0_P = 40.0 / 3.6;        //[m/s]
//		tra.v1_P = 60.0 / 3.6;		//[m/s]
//		tra.v2_P = 60.0 / 3.6;		//[m/s]
//		tra.v0_R = 52.0 / 3.6;        //[m/s]
//		tra.v1_R = 96.0 / 3.6;		//[m/s]
//		tra.v2_R = 96.0 / 3.6;		//[m/s]
//	}
//
//	tra.Fmax = 18.818 * 1000.0;         //[N/MM]  �����p(100%��Ԏ�)(�Q�l�F�N�������x0.639[m/s/s])
//	tra.Bmax = 15.556 * 1000.0;		//[N/MM] (100%��Ԏ�)
//	tra.Bref = 1.111 * 5.0 / 7.0;
//	tra.BEmax = 0.0;
//
//	tra.i1d_P = 90.42780281;  //�͍s���C��g���N�̈�ł�d���d��
//	tra.i1q_P = 237.1149987; //�͍s���C��g���N�̈�ł�q���d��
//	tra.i1d_R = 91.09231795;  //�񐶎��C��g���N�̈�ł�d���d��
//	tra.i1q_R = -211.1226486; //�񐶎��C��g���N�̈�ł�q���d��
//
//	if (tra.Type == 2) {
//		tra.notch = 4.0;
//		tra.Fmax = 14.776 * 1000.0;         //[N/MM]  �����p(100%��Ԏ�)(�Q�l�F�N�������x0.639[m/s/s])
//		tra.Bmax = 13.910 * 1000.0;		//[N/MM] (100%��Ԏ�)
//		tra.Bref = 1.111 * 5.0 / 7.0;
//		tra.BEmax = 0.0;
//
//		tra.i1d_P = 91.877;  //�͍s���C��g���N�̈�ł�d���d��
//		tra.i1q_P = 183.246; //�͍s���C��g���N�̈�ł�q���d��
//		tra.i1d_R = 91.60374;  //�񐶎��C��g���N�̈�ł�d���d��
//		tra.i1q_R = -173.0142; //�񐶎��C��g���N�̈�ł�q���d��
//	}
//
//	tra.B = 1.111 * 5.0 / 7.0;
//
//	tra.direction = ini->direction;
//	tra.UP_DOWN = (double)tra.direction;
//	tra.flag_station = 0;
//	tra.startflag = 0;
//	tra.lap = 0;
//	tra.lapmax = ini->lapmax;
//
//	tra.t_count = 0;
//	tra.c_count = 0;
//
//	tra.wr = 0.0;
//	tra.ws = 0.0;
//	tra.we = 0.0;
//
//	tra.con_speed = 0.0;
//
//	tra.Gr = 6.31;  //�M�A��
//	tra.Ts = 0.00001;  //��������
//	tra.Tf = 0.0001;  //�����萔
//	tra.POLE = 2.0;  //�ɑΐ�
//	tra.rd = 0.41;   //�ԗ֔��a[m]
//	tra.R1 = 0.0970;
//	tra.R2 = 0.07327;  //2����R[��]
//	tra.L1 = 0.030549;
//	tra.L2 = 0.030549;  //2�����ȃC���_�N�^���X[H]
//	tra.M_IM = 0.029514;
//	tra.g0 = 0.001614;
//
//	tra.Cfc = C_FC * ini->num_M;
//	tra.Lfl = L_FL / ini->num_M;
//	tra.Rfl = R_FL / ini->num_M;
//
//	tra.theta = 1.0;
//	tra.theta_old = 1.0;
//
//
//	tra.Emot_pow = 0.0;			//���[�^�͍s�G�l���M�[[J]
//	tra.Emot_reg = 0.0;			//���[�^�񐶃G�l���M�[[J]
//	tra.Emot_air = 0.0;			//�@�B�u���[�L�G�l���M�[[J]
//
//	tra.Emot_pow_c = 0.0;		//���[�^�͍s�G�l���M�[[kWh]
//	tra.Emot_reg_c = 0.0;		//���[�^�񐶃G�l���M�[[kWh]
//	tra.Emot_air_c = 0.0;		//�@�B�u���[�L�G�l���M�[[kWh]
//
//	tra.Eloss_inv_pow = 0.0;
//	tra.Eloss_mot_pow = 0.0;
//	tra.Eloss_fl_pow = 0.0;
//	tra.Eloss_all_pow = 0.0;
//
//	tra.Eloss_inv_reg = 0.0;
//	tra.Eloss_mot_reg = 0.0;
//	tra.Eloss_fl_reg = 0.0;
//	tra.Eloss_all_reg = 0.0;
//
//	tra.Eres_pow = 0.0;
//	tra.Eres_reg = 0.0;
//	tra.Eres_coa = 0.0;
//
//	tra.Esiv_pow = 0.0;
//	tra.Esiv_reg = 0.0;
//	tra.Esiv_coa = 0.0;
//	tra.Esiv_stp = 0.0;
//
//	tra.Eloss_inv_pow_c = 0.0;
//	tra.Eloss_mot_pow_c = 0.0;
//	tra.Eloss_fl_pow_c = 0.0;
//	tra.Eloss_all_pow_c = 0.0;
//
//	tra.Eloss_inv_reg_c = 0.0;
//	tra.Eloss_mot_reg_c = 0.0;
//	tra.Eloss_fl_reg_c = 0.0;
//	tra.Eloss_all_reg_c = 0.0;
//
//	tra.Eres_pow_c = 0.0;
//	tra.Eres_reg_c = 0.0;
//	tra.Eres_coa_c = 0.0;
//
//	tra.Esiv_pow_c = 0.0;
//	tra.Esiv_reg_c = 0.0;
//	tra.Esiv_coa_c = 0.0;
//	tra.Esiv_stp_c = 0.0;
//
//
//
//	tra.Evh = 0.0;				//��ԓ��o�̓G�l���M�[[J]
//	tra.Evh_in = 0.0;			//��ԓ��̓G�l���M�[[J]
//	tra.Evh_out = 0.0;			//��ԏo�̓G�l���M�[[J]
//	tra.Evh_st = 0.0;
//
//	tra.Evh_c = 0.0;			//��ԓ��o�̓G�l���M�[[kWh]
//	tra.Evh_in_c = 0.0;			//��ԓ��̓G�l���M�[[kWh]
//	tra.Evh_out_c = 0.0;		//��ԏo�̓G�l���M�[[kWh]
//	tra.Evh_st_c = 0.0;
//
//	tra.accelflag = 4;
//	tra.brakeflag = 0;
//	tra.laststopflag = 0;
//	tra.X_recentstop = 0;
//	tra.X_nextstop = 0.0;
//	tra.X_brake = 0.0;
//	tra.x = ini->x;
//	tra.v = ini->v;
//	tra.v_new = 0.0;
//	tra.a = 0.0;
//	tra.Fmot = 0.0;
//	tra.Ftot = 0.0;
//
//	tra.id_mot_P = 0.0;
//	tra.id_mot_R = 0.0;
//	tra.iq_mot_P = 0.0;
//	tra.iq_mot_R = 0.0;
//
//	tra.v1d = 0.0;
//	tra.v1q = 0.0;
//	tra.Em = 0.0;
//	tra.idtot = 0.0;
//	tra.iqtot = 0.0;
//	tra.idtotF = 0.0;
//	tra.iqtotF = 0.0;
//	tra.iq_reg = 0.0;
//	tra.iq_regF = 0.0;
//
//	tra.x_c = ini->x / 1000.0;
//	tra.v_c = ini->v / 1000.0;
//
//	tra.Pmot = 0.0;
//	tra.Ploss_mot = 0.0;
//	tra.Ploss_fe = 0.0;
//	tra.Ploss_inv = 0.0;
//	tra.Ploss_fl = 0.0;
//	tra.iinv = 0.0;
//	tra.ifc = 0.0;
//	tra.ifl = 0.0;
//	tra.vfc_old = Vss1;
//	tra.isiv = 0.0;
//	tra.vfc = Vss1;
//	tra.vp = Vss1;
//	tra.vp_old = Vss1;
//	tra.R_grade = 0.0;
//	tra.R_curve = 0.0;
//	tra.R_run = 0.0;
//	tra.R_total = 0.0;
//	tra.t_stop = 0.0;
//	tra.t_stop_old = 0.0;
//	tra.t_wait = 0.0;
//	tra.diaflag = 0;
//
//	tra.Speedlimit = 0.0;
//	tra.Reaccelspeed = 0.0;
//	tra.T_delay = ini->T_delay;
//	tra.Vst = Vss1;
//
//	for (m = 0; m < rownum; m++)
//	{
//		for (n = 0; n < 2; n++)
//		{
//			tra.route[m][n] = 0;
//		}
//	}
//}
//
//void Make_substation(SUB& sub, const INI_SUB* ini, int i)
//{
//
//	//*sub = (SUB*)malloc(sizeof(SUB));
//	//if (*sub == NULL)
//	//{
//	//	//�G���[����
//	//}
//
//	sub.name_SS = i;
//	sub.Vss_0 = Vss1;
//	sub.Iss_0 = Iss_rated1;
//	sub.e_ss = ini->e_ss;
//	sub.Xss = ini->Xss;
//
//	sub.diode = ini->diode;
//	sub.vss = ini->vss;
//	sub.vss_e = ini->vss;
//	sub.iss = ini->iss;
//	sub.vout = ini->vss;
//	sub.flag = 1;
//
//	sub.Rss = (ini->e_ss * ini->vss) / (100.0 * ini->Iss);
//	sub.Rss_ini = (ini->e_ss * ini->vss) / (100.0 * ini->Iss);
//	sub.Jss_ini = -ini->vss / ((ini->e_ss * ini->vss) / (100.0 * ini->Iss));	//Jss�̏����l
//	sub.Jss = -ini->vss / ((ini->e_ss * ini->vss) / (100.0 * ini->Iss));		//�X�V�p��Jss
//
//	sub.Ess = 0.0;			//�ϓd���o�̓G�l���M�[[J]
//	sub.Ess_c = 0.0;		//�ϓd���o�̓G�l���M�[[kWh]
//	sub.Wss = 0.0;
//	sub.Wss_c = 0.0;
//}
//
//void Make_NODE_TRAIN(NODE& node, TRAIN& tra)
//{
//	//*node = (NODE*)malloc(sizeof(NODE));
//	//if (*node == NULL)
//	//{
//	//	printf("MEMORY?? What is that????\n");
//	//	exit(EXIT_FAILURE);//�G���[����
//	//}
//
//	node.Number = tra.name_train;
//	node.X = tra.x;
//	node.V = tra.vp;
//	node.I = tra.ifl;
//	node.r = 0.0;
//	node.flag = tra.direction;
//}
//
//void Make_NODE_SS(NODE& node, SUB& sub)
//{
//	//*node = (NODE*)malloc(sizeof(NODE));
//	//if (*node == NULL)
//	//{
//	//	printf("MEMORY?? What is that????\n");
//	//	exit(EXIT_FAILURE);//�G���[����
//	//}
//
//	node.Number = sub.name_SS;
//	node.X = sub.Xss;
//	node.V = sub.vss;
//	node.I = sub.Jss;
//	node.r = sub.Rss;
//	node.flag = 0;
//}
//
//void Error_Detection_SS(SUB& sub)
//{
//	if (sub.vout > sub.vss_e + 0.1) {
//		sub.Jss = sub.iss - sub.vout / sub.Rss;
//	}
//	else if (sub.vout < sub.vss_e - 0.1) sub.Jss = sub.Jss_ini;
//
//	if (sub.iss < -0.1)		//�ϓd���ɓd�����������Ă����flag=0�i�G���[�j
//	{
//		sub.flag = 0;
//		sub.Jss = -sub.vout / sub.Rss;
//	}
//	else
//	{
//		sub.flag = 1;			//�ϓd���ɓd�����������Ă��Ȃ���flag=1�iOK�j
//	}
//
//	/*	if (sub.vout >(Vcmax))	//�d���̐���
//		{
//			sub.flag = 0;
//			sub.Jss = -(Vcmax) / sub.Rss;
//		}*/
//}
//
//void Error_Detection_REG(TRAIN& tra) {
//
//	/********** �y���׉񐶐���p�^�[���𖞂������ǂ����̔��� **********/
//	if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
//	{
//		tra.b_reg_judge = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
//	}
//	else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
//	{
//		tra.b_reg_judge = 0.0;
//	}
//	else															//�񐶍i�荞�݂Ȃ�
//	{
//		tra.b_reg_judge = tra.Fmot;
//	}
//
//	tra.judgement = fabs(tra.b_reg - tra.b_reg_judge);
//
//	if (tra.judgement > 1.0) {
//		tra.reg_flag = 0;
//	}
//	else {
//		tra.reg_flag = 1;
//	}
//}
//
//void Calculate_R_run(TRAIN& tra)		//���s��R�̎Z�o
//{
//	tra.R_run = (tra.mass_M / 1000.0 * (1.65 + 0.0247 * tra.v_c) + tra.mass_T / 1000.0 * (0.78 + 0.0028 * tra.v_c) + (0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c) * g;		//rr=g�~W�~(a+bv+cv^2/W)
//	//	tra.R_run = g*(1.32+0.0614 * tra.v_c)*tra.mass_c+(0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c;
//}
//
//void Calculate_R_total(TRAIN& tra)		//���v��R�̎Z�o
//{
//	Calculate_R_run(tra);
//
//	tra.R_total = tra.R_run;
//}
//
//
//void Traction_force(TRAIN& tra)		//�͍s�������͂̎Z�o
//{
//
//	double v1_pp = tra.v1_P * tra.vfc / 1350.0;
//	double v2_pp = tra.v2_P * tra.vfc / 1350.0;
//
//	if (tra.Type == 2) tra.Fmax = 14.776 * 1000.0 * (tra.notch / 4.0);
//
//	if (tra.v < v1_pp)		//��g���N�̈�
//	{
//		tra.Fmot = tra.Fmax * tra.num_mot;
//	}
//	else if (tra.v < v2_pp)		//��d�͗̈�
//	{
//		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp / tra.v;
//	}
//	else     //�����̈�
//	{
//		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp * v2_pp / (tra.v * tra.v);	//�����̈�
//	}
//
//	if (tra.Fmot > tra.Fmax * tra.num_mot) tra.Fmot = tra.Fmax * tra.num_mot;
//
//	tra.Ftot = tra.Fmot;
//}
//
//void Regenerative_force(TRAIN& tra)		//�񐶎������͂̎Z�o
//{
//	double v1_rr = tra.v1_R * tra.vfc / 1650.0;
//	double v2_rr = tra.v2_R * tra.vfc / 1650.0;
//
//	if (tra.Type == 2) {
//		tra.Bmax = 13.910 * 1000.0 * tra.brake / 5.0;
//		if (tra.brake > 5.0) {
//			tra.Bmax = 15.092 * 1000.0;
//		}
//	}
//
//	if (tra.v < v1_rr)		//��g���N�̈�
//	{
//		tra.Fmot_regM = -tra.Bmax * tra.num_mot;
//	}
//	else if (tra.v < v2_rr)		//��d�͗̈�
//	{
//		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr / tra.v;
//	}
//	else	//�����̈�
//	{
//		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr * v2_rr / (tra.v * tra.v);
//	}
//
//	tra.Fmot = tra.Fmot_regM;
//	if (tra.Fmot < -tra.BEmax) tra.Fmot = -tra.BEmax;
//	tra.Ftot = -tra.BEmax;
//}
//
///*** �ꎟ�x��t�B���^ ***/
//double funcdelay(double aTs, double aTf, double ayold, double au)
//{
//	double fudely;
//	fudely = (aTs * au + aTf * ayold) / (aTf + aTs);
//	return(fudely);
//}
///************************/
//
//void d_current_power(TRAIN& tra)//�͍s����d���d���Z�o
//{
//	double v0_PP = tra.v0_P * tra.vfc / 1350.0;
//	tra.id_mot_P = tra.i1d_P;
//
//	if (tra.v > v0_PP)		//��ߎ����̈�
//	{
//		tra.id_mot_P = tra.i1d_P * v0_PP / tra.v;
//	}
//
//	tra.idtot = tra.id_mot_P;
//}
//
//void d_current_regenerative(TRAIN& tra)//�񐶎���d���d���Z�o
//{
//	double v0_RR = tra.v0_R * tra.vfc / 1650.0;
//	tra.id_mot_R = tra.i1d_R;
//
//	if (tra.v > v0_RR)		//��ߎ����̈�
//	{
//		tra.id_mot_R = tra.i1d_R * v0_RR / tra.v;
//	}
//	tra.idtot = tra.id_mot_R;
//}
//
//void q_current_power(TRAIN& tra)//�͍s����q���d���Z�o
//{
//	tra.iq_mot_P = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
//	tra.iqtot = tra.iq_mot_P;
//	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
//}
//
//void q_current_regenerative(TRAIN& tra)//�񐶎���q���d���Z�o
//{
//	tra.iq_mot_R = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
//	tra.iqtot = tra.iq_mot_R;
//	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
//}
//
//void d_q_voltage_calculation(TRAIN& tra) {
//	if (tra.accelflag == 2) {
//		//�p���g���v�Z
//		tra.ws = 0.0;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = 0.0;
//
//		/*q���d��*/
//		tra.v1q = 0.0;
//		tra.Em = 0.0;
//	}
//	else if (tra.accelflag == 4) {
//
//		tra.ws = 0.0;
//		tra.we = 0.0;
//		/*d���d��*/
//		tra.v1d = 0.0;
//
//		/*q���d��*/
//		tra.v1q = 0.0;
//		tra.Em = 0.0;
//	}
//	else {
//		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
//		//		tra.v1d = - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtot;
//
//				/*q���d��*/
//		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
//		//		tra.v1q = (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
//
//				/*�C���o�[�^���͓d��*/
//		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
//	}
//}
//
//void Solve_next_state(TRAIN& tra)
//{
//	tra.a = (tra.Ftot - tra.R_total) / (tra.mass);		//�^��������ma=F
//	tra.a_c = tra.a * 3600.0 / 1000.0;
//
//	tra.v_new = tra.v + tra.a * dt;
//
//	tra.v = tra.v_new;
//	tra.v_c = tra.v * 3.6;
//
//	tra.wr = (tra.v / tra.rd) * tra.Gr * tra.POLE;
//
//	if (tra.direction == 1)
//	{
//		tra.x += tra.v * dt;
//	}
//	else
//	{
//		tra.x -= tra.v * dt;
//	}
//
//	tra.x_c = tra.x / 1000.0;
//}
//
//void Calculate_BEmax(TRAIN& tra)
//{
//	if (tra.Type == 2) tra.B = 1.111 * tra.brake / 7.0;
//	tra.BEmax = tra.B * tra.mass;
//}
//
//void change_direction(TRAIN& tra)
//{
//	tra.t_count = 0;
//	tra.c_count = 0;
//
//	switch (tra.direction) {
//	case 1:
//		tra.direction = -1;
//		tra.UP_DOWN = -1;
//		break;
//	case -1:
//		tra.direction = 1;
//		tra.UP_DOWN = 1;
//		break;
//	}
//}
//
//void decide_final_station(TRAIN& tra, int i)
//{
//	if (tra.Type == 0 || tra.Type == 9 || tra.Type == 10 || tra.Type == 11 || tra.Type == 12 || tra.Type == 13 || tra.Type == 14) {
//		if (tra.direction == 1.0)			//�r�܁ˏa�J
//		{
//			if (i == NUM_station)
//			{
//				if (tra.lap != tra.lapmax) {
//					change_direction(tra);
//					if (tra.name_train == 7) tra.T_delay = 27.0 * 60.0;
//					if (tra.name_train == 8) tra.T_delay = 37.0 * 60.0;
//					tra.lap++;
//				}
//				else {
//					tra.laststopflag = 1;
//					tra.T_delay = 9999.0;
//				}
//			}
//		}
//		else
//		{
//			if (i == 0)						//�a�J�˒r��
//			{
//				if (tra.lap != tra.lapmax) {
//					change_direction(tra);
//					if (tra.name_train == 3) tra.T_delay = 27.0 * 60.0 + 50.0;
//					if (tra.name_train == 4) tra.T_delay = 37.0 * 60.0 + 30.0;
//					tra.lap++;
//				}
//				else {
//					tra.laststopflag = 1;
//					tra.T_delay = 9999.0;
//				}
//			}
//		}
//	}
//	else if (tra.Type == 1) {
//		if (tra.direction == 1.0)			//��{�ː������
//		{
//			if (i == NUM_station2)
//			{
//				if (tra.lap != tra.lapmax) {
//					change_direction(tra);
//					tra.lap++;
//				}
//				else {
//					tra.laststopflag = 1;
//					tra.T_delay = 9999.0;
//				}
//			}
//		}
//		else
//		{
//			if (i == 0)						//����ˑ�{����
//			{
//				if (tra.lap != tra.lapmax) {
//					change_direction(tra);
//					if (tra.name_train == 5) tra.T_delay = 41.0 * 60.0 + 50.0;
//					tra.lap++;
//				}
//				else {
//					tra.laststopflag = 1;
//					tra.T_delay = 9999.0;
//				}
//			}
//		}
//	}
//}
//
//
////������ς�����by Kasai1012
//void Run_pattern(TRAIN& tra, const STATION* sta, double t, std::vector<double> wait_time)
//{
//	int i;
//
//	tra.t_stop = t;
//
//	/*** ���̉w��T�� ***/
//	tra.X_nextstop_old = tra.X_nextstop;
//	if (tra.Type == 0 || tra.Type == 9 || tra.Type == 10 || tra.Type == 11 || tra.Type == 12 || tra.Type == 13 || tra.Type == 14) {
//		for (i = 0; i < NUM_station + NUM_final_station; i++)
//		{
//			if (tra.direction == 1.0)			//�r�܁ˏa�J
//			{
//				if (sta[i].Xs <= tra.x && sta[i + 1].Xs >= tra.x)
//				{
//					tra.X_nextstop = sta[i + 1].Xs;
//					break;
//				}
//			}
//			else								//�a�J�˒r��
//			{
//				if (sta[i].Xs <= tra.x && sta[i + 1].Xs >= tra.x)
//				{
//					tra.X_nextstop = sta[i].Xs;
//					break;
//				}
//			}
//		}
//
//		decide_final_station(tra, i);
//
//	}
//
//	/*���x�����ݒ�*/
//	if (sta[1].Xs < tra.x && tra.x < sta[2].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[2].Xs < tra.x && tra.x < sta[3].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[3].Xs < tra.x && tra.x < sta[4].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[4].Xs < tra.x && tra.x < sta[5].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[5].Xs < tra.x && tra.x < sta[6].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[6].Xs < tra.x && tra.x < sta[7].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else if (sta[7].Xs < tra.x && tra.x < sta[8].Xs)
//	{
//		tra.Speedlimit = 65.0;
//	}
//	else
//	{
//		tra.Speedlimit = 65.0;
//	}
//
//
//	/**** ���s��R�v�Z ****/
//	Calculate_R_total(tra);
//
//	tra.Bref = 1.111 * tra.brake / 7.0;
//
//
//
//	/*** ���s���[�h������ ***/
//	if (tra.accelflag == 4)		/*��Ԓ��̓��������*/
//	{
//		if (tra.laststopflag == 1)
//		{
//			tra.accelflag = 4;
//		}
//		else if (t < tra.T_delay)
//		{
//			tra.accelflag = 4;
//		}
//		else if (tra.direction == 1) {	//�����i�r�܁��a�J�j
//			if (sta[2].Xs - 200.0 <= tra.x && tra.x < sta[2].Xs + 200.0) {	//�ڔ�
//				// if (tra.t_stop >= tra.T_delay + 2.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;�@�@���Ԏ�����ݒ肷��ꍇ
//				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[0])) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[3].Xs - 200.0 <= tra.x && tra.x < sta[3].Xs + 200.0) {	//���c�n��
//				// if (tra.t_stop >= tra.T_delay + 4.0 * 60.0 + 20 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[0] + wait_time[1])) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[4].Xs - 200.0 <= tra.x && tra.x < sta[4].Xs + 200.0) {	//�V��v��
//				// if (tra.t_stop >= tra.T_delay + 6.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[1] + wait_time[2])) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[5].Xs - 200.0 <= tra.x && tra.x < sta[5].Xs + 200.0) {	//�V�h
//				// if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[6].Xs - 200.0 <= tra.x && tra.x < sta[6].Xs + 200.0) {	//��X��
//				// if (tra.t_stop >= tra.T_delay + 10.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[3])) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[7].Xs - 200.0 <= tra.x && tra.x < sta[7].Xs + 200.0) {	//���h
//				// if (tra.t_stop >= tra.T_delay + 13.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1; 
//				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
//				if (tra.Type == 0) {
//					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[3] + wait_time[4])) tra.accelflag = 1;
//				}
//
//			}
//			else if (sta[8].Xs - 200.0 <= tra.x && tra.x < sta[8].Xs + 200.0) {	//�a�J
//				if ((tra.t_stop - tra.t_stop_old) >= 7200.0) tra.accelflag = 1;
//			}
//			else {
//				if ((tra.t_stop - tra.t_stop_old) >= 20.0) tra.accelflag = 1;
//			}
//		}
//
//
//		else
//		{
//			tra.accelflag = 4;
//		}
//
//	}
//	else if (tra.accelflag == 2)		/*�čs���̓��������*/
//	{
//		if (tra.direction == 1.0)		//��{�ː������
//		{
//			if (tra.x + tra.v * tra.v / 2.0 / tra.Bref < tra.X_nextstop)
//			{
//				tra.accelflag = 2;
//			}
//			/*	else if(tra.v_c < 50.0)
//				{
//					tra.accelflag = 1;
//				}
//				*/
//			else
//			{
//				tra.accelflag = 3;
//				tra.brakeflag = 1;
//			}
//		}
//		else							//��{�ː������
//		{
//			if (tra.x - tra.v * tra.v / 2.0 / tra.Bref > tra.X_nextstop)
//			{
//				tra.accelflag = 2;
//			}
//
//			/*else if(tra.v_c < 60.0)
//			{
//				tra.accelflag = 1;
//			}
//			*/
//			else
//			{
//				tra.accelflag = 3;
//				tra.brakeflag = 1;
//			}
//		}
//	}
//	else		/*������or�������̓��������*/
//	{
//		if (tra.direction == 1.0)								//��{�ː������
//		{
//			if (tra.brakeflag == 1)		/*�������̓���*/
//			{
//				if (tra.v_c >= 0.0)
//				{
//					tra.accelflag = 3;		//�������
//				}
//				else
//				{
//					tra.accelflag = 4;		//���
//					tra.brakeflag = 0;
//					tra.t_stop_old = tra.t_stop;
//				}
//			}
//			else
//			{
//				if (tra.laststopflag == 1)	//�w�̎�O�ɒ�Ԃ����ꍇ�̋~��
//				{
//					tra.accelflag = 3;
//					tra.brakeflag = 1;
//				}
//				else if (tra.x + tra.v * tra.v / 2.0 / tra.Bref > tra.X_nextstop && tra.v_c > speed_min - 20.0)
//				{
//					tra.accelflag = 3;
//				}
//				else if ((tra.v_c <= tra.Speedlimit) && (tra.a >= 0))
//				{
//					tra.accelflag = 1;		//�������
//				}
//				else
//				{
//					tra.accelflag = 2;		//�čs���
//				}
//			}
//		}
//		else													//����ˑ�{����
//		{
//			if (tra.brakeflag == 1)		/*�������̓���*/
//			{
//				if (tra.v_c >= 0.0)
//				{
//					tra.accelflag = 3;		//�������
//				}
//				else
//				{
//					tra.accelflag = 4;		//���
//					tra.brakeflag = 0;
//					tra.t_stop_old = tra.t_stop;
//				}
//			}
//			else
//			{
//				if (tra.laststopflag == 1)	//�w�̎�O�ɒ�Ԃ����ꍇ�̋~��
//				{
//					tra.accelflag = 3;
//					tra.brakeflag = 1;
//				}
//				else if (tra.x - tra.v * tra.v / 2.0 / tra.Bref < tra.X_nextstop && tra.v_c > speed_min - 20.0)
//				{
//					tra.accelflag = 3;
//				}
//				else if ((tra.v_c <= tra.Speedlimit) && (tra.a >= 0))
//				{
//					tra.accelflag = 1;		//�������
//				}
//				else
//				{
//					tra.accelflag = 2;		//�čs���
//				}
//			}
//		}
//	}
//
//}
//
//void Calculation_traction_force_brake_force(TRAIN& tra) {
//
//	/*** �e���s���[�h�ɉ�����������w�� ***/
//	if (tra.accelflag == 1)				//������
//	{
//		/** �u���[�L�͌v�Z(�������Ȃ̂�0) **/
//		tra.b_reg = 0.0;
//		tra.b_air = 0.0;
//		tra.b_reg_loss = 0.0;
//
//		/** �ԗ������͌v�Z **/
//		Traction_force(tra);
//
//		/**�͍s��d���d���v�Z**/
//		d_current_power(tra);
//
//		/**�͍s��q���d���v�Z**/
//		q_current_power(tra);
//
//		tra.iq_reg = 0.0;
//
//		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
//
//		/*q���d��*/
//		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
//		/*�C���o�[�^���͓d��*/
//		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
//	}
//	else if (tra.accelflag == 2)			//�čs��
//	{
//		/** �u���[�L�͌v�Z(�čs���Ȃ̂�0) **/
//		tra.b_reg = 0.0;
//		tra.b_air = 0.0;
//		tra.b_reg_loss = 0.0;
//
//		/** �ԗ������͌v�Z **/
//		tra.Fmot = 0.0;
//		tra.Ftot = 0.0;
//
//		/*�d���v�Z*/
//		tra.idtot = 0.0;
//		tra.iqtot = 0.0;
//		tra.iqtotF = 0.0;
//		tra.iq_reg = 0.0;
//
//		//�p���g���v�Z
//		tra.ws = 0.0;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = 0.0;
//
//		/*q���d��*/
//		tra.v1q = 0.0;
//		tra.Em = 0.0;
//
//	}
//	else if (tra.accelflag == 3)			//������
//	{
//		/** �ԗ��u���[�L�͌v�Z **/
//		Regenerative_force(tra);
//		d_current_regenerative(tra);
//		q_current_regenerative(tra);
//
//		/********** �y���׉񐶒��̃g���N�i�荞�ݖ͋[���i�y���׉񐶐��䉺�̓d�C�u���[�L�͂��v�Z�j **********/
//		if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
//		{
//			tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
//			tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
//			tra.iqtot = tra.iq_reg;
//		}
//		else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
//		{
//			tra.b_reg = 0.0;
//			tra.iq_reg = 0.0;
//			tra.iqtot = 0.0;
//		}
//		else															//�񐶍i�荞�݂Ȃ�
//		{
//			tra.b_reg = tra.Fmot;
//			tra.iq_reg = tra.iq_mot_R;
//			tra.iqtot = tra.iq_reg;
//		}
//
//		tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
//		tra.iqtotF = tra.iq_regF;
//
//		tra.b_reg_loss = tra.Fmot - tra.b_reg;		//�񐶍i�荞�ݗ�[N]
//
//		tra.b_air = -tra.BEmax - tra.b_reg;			//����Ȃ��u���[�L�͂͋�C�u���[�L�ŕ₤
//
//		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
//
//		/*q���d��*/
//		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
//		tra.Em = sqrt(tra.v1q * tra.v1q + tra.v1d * tra.v1d);
//	}
//	else if (tra.accelflag == 4)						//��Ԓ�
//	{
//		tra.Fmot = 0.0;
//		tra.Ftot = 0.0;
//		tra.idtot = 0.0;
//		tra.iqtot = 0.0;
//		tra.iqtotF = 0.0;
//		tra.v = 0.0;
//		tra.v_c = 0.0;
//		tra.a = 0.0;
//		tra.b_reg = 0.0;
//		tra.b_air = 0.0;
//		tra.b_reg_loss = 0.0;
//		tra.iq_reg = 0.0;
//
//		tra.ws = 0.0;
//		tra.we = 0.0;
//		tra.wr = 0.0;
//		/*d���d��*/
//		tra.v1d = 0.0;
//
//		/*q���d��*/
//		tra.v1q = 0.0;
//		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
//	}
//	else if (tra.accelflag == 5) {
//		if (tra.v_c > tra.con_speed) {
//			tra.Fmot = -tra.Bmax * tra.num_mot * (tra.v_c - tra.con_speed) / 5.0;
//			tra.Ftot = tra.Fmot;
//
//			d_current_regenerative(tra);
//			q_current_regenerative(tra);
//
//			/********** �y���׉񐶒��̃g���N�i�荞�ݖ͋[���i�y���׉񐶐��䉺�̓d�C�u���[�L�͂��v�Z�j **********/
//			if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
//			{
//				tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
//				tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
//				tra.iqtot = tra.iq_reg;
//			}
//			else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
//			{
//				tra.b_reg = 0.0;
//				tra.iq_reg = 0.0;
//				tra.iqtot = 0.0;
//			}
//			else															//�񐶍i�荞�݂Ȃ�
//			{
//				tra.b_reg = tra.Fmot;
//				tra.iq_reg = tra.iq_mot_R;
//				tra.iqtot = tra.iq_reg;
//			}
//
//			tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
//			tra.iqtotF = tra.iq_regF;
//
//			tra.b_reg_loss = tra.Fmot - tra.b_reg;		//�񐶍i�荞�ݗ�[N]
//
//			tra.b_air = tra.Ftot - tra.b_reg;			//����Ȃ��u���[�L�͂͋�C�u���[�L�ŕ₤
//		}
//		else {
//			tra.Fmot = tra.Fmax * tra.num_mot * (tra.con_speed - tra.v_c) / 5.0;
//			tra.Ftot = tra.Fmot;
//
//			d_current_power(tra);
//			q_current_power(tra);
//
//			tra.iq_reg = 0.0;
//		}
//
//		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
//		tra.we = tra.ws + tra.wr;
//
//		/*d���d��*/
//		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
//
//		/*q���d��*/
//		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
//		tra.Em = sqrt(tra.v1q * tra.v1q + tra.v1d * tra.v1d);
//	}
//
//	tra.Fmot_c = tra.Fmot / 1000.0;
//	tra.b_reg_c = tra.b_reg / 1000.0;
//	tra.Ftot_c = tra.Ftot / 1000.0;
//	tra.b_air_c = tra.b_air / 1000.0;
//	tra.b_reg_loss_c = tra.b_reg_loss / 1000.0;
//}
//
//void Calculation_traction_circuit(TRAIN& tra)
//{
//	tra.isiv = P_SIV * tra.num / tra.vfc;
//
//	if (tra.Type == 3) tra.isiv = P_KARAKIDA / tra.vfc;
//	if (tra.Type == 4) tra.isiv = P_SHINYURIGAOKA / 1500.0;
//
//	tra.ifl = tra.isiv + (tra.Pmot + tra.Ploss_mot + tra.Ploss_fe) / tra.vfc;
//
//	tra.Ploss_mot = (tra.R1 * (tra.idtot * tra.idtot + tra.iqtotF * tra.iqtotF) + tra.R2 * (tra.M_IM / tra.L2 * tra.iqtotF) * (tra.M_IM / tra.L2 * tra.iqtotF)) * tra.num_mot;//���[�^�̓���
//	tra.Ploss_fe = tra.g0 * (tra.v1d * tra.v1d + tra.v1q * tra.v1q) * tra.num_mot;
//	tra.Ploss_fl = tra.Rfl * tra.ifl * tra.ifl;					//�t�B���^���A�N�g���̃W���[������[W]
//	tra.Ploss_fl_c = tra.Ploss_fl / 1000.0;					//�t�B���^���A�N�g���̃W���[������[kW]
//
//	tra.Ploss_mot_c = (tra.Ploss_mot + tra.Ploss_fe) / 1000.0;				//���[�^�̓���[kW]
//
//	tra.Ploss_all = tra.Ploss_mot + tra.Ploss_fe + tra.Ploss_fl;
//	tra.Ploss_all_c = tra.Ploss_all / 1000.0;
//}
//
//void Calculation_power_motor(TRAIN& tra)
//{
//	tra.Ptot = tra.v * tra.Ftot;
//
//	if (tra.accelflag == 1)				/*�͍s��*/
//	{
//		tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
//		tra.Pmot_air = 0.0;
//		tra.Pmot_reg_loss = 0.0;
//
//	}
//	else if (tra.accelflag == 3)		/*�񐶒�*/
//	{
//		tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
//		tra.Pmot_air = tra.v * tra.b_air;
//		tra.Pmot_reg_loss = tra.v * tra.b_reg_loss;
//	}
//	else if (tra.accelflag == 5) {
//		if (tra.v_c > tra.con_speed) {
//			tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
//			tra.Pmot_air = tra.v * tra.b_air;
//			tra.Pmot_reg_loss = tra.v * tra.b_reg_loss;
//		}
//		else {
//			tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
//			tra.Pmot_air = 0.0;
//			tra.Pmot_reg_loss = 0.0;
//		}
//	}
//	else
//	{
//		tra.Pmot = 0.0;
//		tra.Pmot_air = 0.0;
//		tra.Pmot_reg_loss = 0.0;
//	}
//
//	tra.Pres = tra.v * tra.R_total;		//���s��R�ł̃��X
//
//	tra.Ptot_c = tra.Ptot / 1000.0;
//	tra.Pmot_c = tra.Pmot / 1000.0;
//	tra.Pmot_air_c = tra.Pmot_air / 1000.0;
//	tra.Pres_c = tra.Pres / 1000.0;
//	tra.Pmot_reg_loss_c = tra.Pmot_reg_loss / 1000.0;
//
//}
//
//void Calculation_power_train(TRAIN& tra)
//{
//	tra.Pvh = tra.vp * tra.ifl;
//
//	if (tra.accelflag == 1 || tra.accelflag == 2)				//�͍s��,�čs��
//	{
//		tra.Pvh_in = tra.vp * tra.ifl;
//		tra.Pvh_out = 0.0;
//		tra.Pvh_st = 0.0;
//	}
//	else if (tra.accelflag == 3)		//�񐶒�
//	{
//		tra.Pvh_in = 0.0;
//		tra.Pvh_out = tra.vp * tra.ifl;
//		tra.Pvh_st = 0.0;
//	}
//	else if (tra.accelflag == 5) {
//		if (tra.v_c > tra.con_speed) {
//			tra.Pvh_in = 0.0;
//			tra.Pvh_out = tra.vp * tra.ifl;
//			tra.Pvh_st = 0.0;
//		}
//		else {
//			tra.Pvh_in = tra.vp * tra.ifl;
//			tra.Pvh_out = 0.0;
//			tra.Pvh_st = 0.0;
//		}
//	}
//	else
//	{
//		tra.Pvh_in = 0.0;
//		tra.Pvh_out = 0.0;
//		tra.Pvh_st = tra.vp * tra.ifl;
//	}
//	tra.Pvh_c = tra.Pvh / 1000.0;
//	tra.Pvh_in_c = tra.Pvh_in / 1000.0;
//	tra.Pvh_out_c = tra.Pvh_out / 1000.0;
//	tra.Pvh_st_c = tra.Pvh_st / 1000;
//}
//
//void Calculation_energy_train(TRAIN& tra)
//{
//	if (tra.accelflag == 1)				//�͍s��
//	{
//		tra.Emot_pow += tra.Pmot * dt;
//
//		tra.Eloss_inv_pow += tra.Ploss_inv * dt;
//		tra.Eloss_mot_pow += (tra.Ploss_mot + tra.Ploss_fe) * dt;
//		tra.Eloss_fl_pow += tra.Ploss_fl * dt;
//		tra.Eloss_all_pow += tra.Ploss_all * dt;
//		tra.Esiv_pow += P_SIV * dt;
//
//		tra.Eres_pow += tra.Pres * dt;
//		//		tra.Emot_pow += (tra.Pvh_in - tra.Ploss_all - P_SIV - tra.Pres)* dt;
//	}
//	else if (tra.accelflag == 3)		//�񐶒�
//	{
//		tra.Emot_reg += tra.Pmot * dt;
//		tra.Emot_reg_loss += tra.Pmot_reg_loss * dt;
//
//		tra.Eloss_inv_reg += tra.Ploss_inv * dt;
//		tra.Eloss_mot_reg += (tra.Ploss_mot + tra.Ploss_fe) * dt;
//		tra.Eloss_fl_reg += tra.Ploss_fl * dt;
//		tra.Eloss_all_reg += tra.Ploss_all * dt;
//		tra.Esiv_reg += P_SIV * dt;
//
//		tra.Eres_reg += tra.Pres * dt;
//		//		tra.Emot_reg += (tra.Pvh_out - tra.Ploss_all - P_SIV - tra.Pres)* dt;
//	}
//	else if (tra.accelflag == 2)
//	{
//		tra.Eres_coa += tra.Pres * dt;
//		tra.Esiv_coa += P_SIV * dt;
//	}
//	else if (tra.accelflag == 4)
//	{
//		tra.Esiv_stp += P_SIV * dt;
//		if (tra.Type == 3) tra.Esiv_stp += P_KARAKIDA * dt;
//		if (tra.Type == 4) tra.Esiv_stp += P_SHINYURIGAOKA * dt;
//	}
//	else if (tra.accelflag == 5) {
//		if (tra.v_c > tra.con_speed) {
//			tra.Emot_reg += tra.Pmot * dt;
//			tra.Emot_reg_loss += tra.Pmot_reg_loss * dt;
//
//			tra.Eloss_inv_reg += tra.Ploss_inv * dt;
//			tra.Eloss_mot_reg += (tra.Ploss_mot + tra.Ploss_fe) * dt;
//			tra.Eloss_fl_reg += tra.Ploss_fl * dt;
//			tra.Eloss_all_reg += tra.Ploss_all * dt;
//			tra.Esiv_reg += P_SIV * dt;
//
//			tra.Eres_reg += tra.Pres * dt;
//		}
//		else {
//			tra.Emot_pow += tra.Pmot * dt;
//
//			tra.Eloss_inv_pow += tra.Ploss_inv * dt;
//			tra.Eloss_mot_pow += (tra.Ploss_mot + tra.Ploss_fe) * dt;
//			tra.Eloss_fl_pow += tra.Ploss_fl * dt;
//			tra.Eloss_all_pow += tra.Ploss_all * dt;
//			tra.Esiv_pow += P_SIV * dt;
//
//			tra.Eres_pow += tra.Pres * dt;
//		}
//	}
//
//	tra.Emot_air += tra.Pmot_air * dt;
//
//
//	tra.Emot_pow_c = tra.Emot_pow / 1000.0 / 3600.0;
//	tra.Emot_reg_c = tra.Emot_reg / 1000.0 / 3600.0;
//	tra.Emot_air_c = tra.Emot_air / 1000.0 / 3600.0;
//	tra.Emot_reg_loss_c = tra.Emot_reg_loss / 1000.0 / 3600.0;
//
//	tra.Eloss_inv_pow_c = tra.Eloss_inv_pow / 1000.0 / 3600.0;
//	tra.Eloss_mot_pow_c = tra.Eloss_mot_pow / 1000.0 / 3600.0;
//	tra.Eloss_fl_pow_c = tra.Eloss_fl_pow / 1000.0 / 3600.0;
//	tra.Eloss_all_pow_c = tra.Eloss_all_pow / 1000.0 / 3600.0;
//
//	tra.Eloss_inv_reg_c = tra.Eloss_inv_reg / 1000.0 / 3600.0;
//	tra.Eloss_mot_reg_c = tra.Eloss_mot_reg / 1000.0 / 3600.0;
//	tra.Eloss_fl_reg_c = tra.Eloss_fl_reg / 1000.0 / 3600.0;
//	tra.Eloss_all_reg_c = tra.Eloss_all_reg / 1000.0 / 3600.0;
//
//	tra.Eres_pow_c = tra.Eres_pow / 1000.0 / 3600.0;
//	tra.Eres_reg_c = tra.Eres_reg / 1000.0 / 3600.0;
//	tra.Eres_coa_c = tra.Eres_coa / 1000.0 / 3600.0;
//
//	tra.Esiv_pow_c = tra.Esiv_pow / 1000.0 / 3600.0;
//	tra.Esiv_reg_c = tra.Esiv_reg / 1000.0 / 3600.0;
//	tra.Esiv_coa_c = tra.Esiv_coa / 1000.0 / 3600.0;
//	tra.Esiv_stp_c = tra.Esiv_stp / 1000.0 / 3600.0;
//
//	tra.Evh += tra.Pvh * dt;
//	tra.Evh_in += tra.Pvh_in * dt;
//	tra.Evh_out += tra.Pvh_out * dt;
//	tra.Evh_st += tra.Pvh_st * dt;
//
//	tra.Evh_c = tra.Evh / 1000.0 / 3600.0;
//	tra.Evh_in_c = tra.Evh_in / 1000.0 / 3600.0;
//	tra.Evh_out_c = tra.Evh_out / 1000.0 / 3600.0;
//	tra.Evh_st_c = tra.Evh_st / 1000.0 / 3600.0;
//}
//
//void Calculation_power_sub(SUB& sub)
//{
//	sub.Pss = sub.iss * sub.vout;
//	sub.Pss_c = sub.Pss / 1000.0;
//}
//
//void Calculation_energy_sub(SUB& sub)
//{
//	if (t >= 15.0 * 60.0 && t < 45.0 * 60.0) {
//		sub.Ess += sub.Pss * dt;						//[J]
//		sub.Ess_c = sub.Ess / 1000.0 / 3600.0;		//[kWh]
//		if (sub.diode == 1)
//		{
//			sub.Wss += sub.Rss * sub.iss * sub.iss * dt;		//[J]
//			sub.Wss_c = sub.Wss / 1000 / 3600;				//[kWh]
//		}
//	}
//}
//
///*** x�����y���w���v�f������ ***/
//void sort_s(const std::unique_ptr<NODE[]>& data, int n)
//{
//
//	std::vector<NODE> temp(N_node);
//	for (size_t i = 0; i < n; i++)temp.at(i) = data[i];
//
//	std::sort(temp.begin(), temp.end(), [](auto const& ll, auto const& rr) {return ll.X < rr.X; });
//
//	for (size_t i = 0; i < n; i++)data[i] = temp[i];
//
//	/*int k = n - 1;
//	while (k >= 0)
//	{
//		int i, j;
//		for (i = 1, j = -1; i <= k; i++)
//			if (data[i - 1].X > data[i].X) {
//				j = i - 1;
//				swap_node(data[i], data[j]);
//			}
//		k = j;
//	}*/
//
//	return;
//}
//
//
///*** �z��data[]�̐擪n�̗v�f�������̏����Ƀ\�[�g ***/
//void sort(const std::unique_ptr<NODE[]>& data, int n)
//{
//	int k = n - 1;
//	while (k >= 0)
//	{
//		size_t i = 0, j = 0;
//		for (i = 1, j = -1; i <= k; i++)
//			if (data[i - 1].X > data[i].X) {
//				j = i - 1;
//				NODE tempI, tempJ;
//				tempI = data[i];
//				tempJ = data[j];
//				data[j] = tempI;
//				data[i] = tempJ;
//				//swap(data[i], data[j]);
//			}
//		k = j;
//	}
//}
//
//void Initialize_BRANCH(BRANCH& branch)
//{
//	//*branch = (BRANCH*)malloc(sizeof(BRANCH));
//	//if (*branch == NULL)
//	//{
//	//	//�G���[����
//	//}
//}
//
//int Count_Trains_Direction(const std::unique_ptr<TRAIN[]>& tra, int direction) {
//	int count_d = 0;
//	for (int i = 0; i < NUM_tra; i++) {
//		if (tra[i].direction == direction) count_d++;
//	}
//	return count_d;
//}
//
//void Make_branch(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<NODE[]>& temp, const std::unique_ptr<TRAIN[]>& tra)
//{
//	int i, j;						//j�F�������ɐڑ�����Ă���m�[�h�T���p�Ck�F�������ɐڑ�����Ă���m�[�h�T���p
//	int count_SS, count_UP, count_DOWN, NUM_tra_UP;	//�u�����`�쐬�p�J�E���g�ϐ�
//
//	count_SS = 0;
//	count_UP = NUM_sub;
//	NUM_tra_UP = Count_Trains_Direction(tra, 1);
//
//	//printf("%d\n", NUM_tra_UP);
//
//	count_DOWN = 2 * NUM_sub + NUM_tra_UP - 1;
//
//	/*�@SS�̃u�����`���쐬*/
//	//printf("-----------SS--------------\n");
//	for (i = 0; i < N_node; i++)
//	{
//		//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
//		if (temp[i].flag == 0)				//In case of SS( i Starts with 0 to NUM_sub-1 )
//		{
//			//puts("SS");
//			branch[count_SS].X_start = temp[i].X;		//�ʒu[m]�i��_����̋����j
//			branch[count_SS].X_end = temp[i].X;		//�ʒu[m]�i��_����̋����j
//			branch[count_SS].r = temp[i].r;			//Resistance [Ohm]
//			branch[count_SS].flag = temp[i].flag;		//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
//			branch[count_SS].Node_pos = -1;			//���������ɐڑ�����Ă���m�[�h��
//			branch[count_SS].Node_neg = -1;			//���������ɐڑ�����Ă���m�[�h��
//			count_SS = count_SS + 1;
//		}
//	}
//
//	/*�A���̃u�����`���쐬*/
//	//printf("-----------UP--------------\n");
//	for (i = 0; i < N_node; i++)
//	{
//		if (temp[i].flag >= 0)			//In case of +1 direction
//		{
//			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
//			if (count_UP < NUM_tra_UP + 2 * NUM_sub - 1)		//�I�[�ȊO
//			{
//
//				for (j = 0; temp[i + (j + 1)].flag == -1; j++) {}				//temp[i]�ɑ΂��Đ��������ɐڑ�����Ă���m�[�h��T��
//				//printf("\nj = %d\n", j);
//
//				/*temp[i]�ɑ΂��Đ������ɐڑ�����Ă���u�����`B���쐬*/
//				branch[count_UP].Node_pos = temp[i + (j + 1)].Number;			//�u�����`B�ɂ��āC���������ɐڑ�����Ă���m�[�h��
//				branch[count_UP].Node_neg = temp[i].Number;					//�u�����`B�ɂ��āC���������ɐڑ�����Ă���m�[�h��
//
//				branch[count_UP].X_start = temp[i].X;							//[m]
//				branch[count_UP].X_end = temp[i + (j + 1)].X;					//[m]
//
//				if (temp[i + (j + 1)].X <= 7500.0) {
//					branch[count_UP].r = R_StoN_UP * fabs(temp[i + (j + 1)].X - temp[i].X) / 1000.0;		//Resistance [Ohm]
//				}
//				else {
//					branch[count_UP].r = R_NtoK_UP * fabs(temp[i + (j + 1)].X - temp[i].X) / 1000.0;		//Resistance [Ohm]
//				}
//
//				branch[count_UP].flag = 1;										//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
//
//				/*���̃��[�v�̂��߂̏���*/
//				j = 0;
//				count_UP = count_UP + 1;
//			}
//			else {}		//�I�[�̃m�[�h�ł͉������Ȃ�
//		}
//	}
//
//
//	/*�B����̃u�����`���쐬*/
//	//printf("-----------DOWN--------------\n");
//	for (i = 0; i < N_node; i++)
//	{
//		if (temp[i].flag <= 0)			//In case of +1 direction
//		{
//			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
//			if (count_DOWN < N_branch)		//�I�[�ȊO
//			{
//				for (j = 0; temp[i + (j + 1)].flag == 1; j++) {}				//temp[i]�ɑ΂��Đ��������ɐڑ�����Ă���m�[�h��T��
//				//printf("\nj = %d\n", j);
//
//				/*�������ɐڑ�����Ă���u�����`���쐬*/
//				branch[count_DOWN].Node_pos = temp[i + (j + 1)].Number;		//���������ɐڑ�����Ă���m�[�h��
//				branch[count_DOWN].Node_neg = temp[i].Number;					//���������ɐڑ�����Ă���m�[�h��
//
//				branch[count_DOWN].X_start = temp[i].X;						//[m]
//				branch[count_DOWN].X_end = temp[i + (j + 1)].X;				//[m]
//
//				if (temp[i + (j + 1)].X <= 7500.0) {
//					branch[count_DOWN].r = R_StoN_DOWN * fabs(temp[i + (j + 1)].X - temp[i].X) / 1000.0;		//Resistance [Ohm]
//				}
//				else {
//					branch[count_DOWN].r = R_NtoK_DOWN * fabs(temp[i + (j + 1)].X - temp[i].X) / 1000.0;		//Resistance [Ohm]
//				}
//
//				branch[count_DOWN].flag = -1;									//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
//
//				/*���̃��[�v�̂��߂̏���*/
//				j = 0;
//				count_DOWN = count_DOWN + 1;
//			}
//			else {}		//�I�[�̃m�[�h�ł͉������Ȃ�
//		}
//	}
//}
//
//void Initialize_matrix1(const std::unique_ptr<double[]>& Mtr, size_t M, size_t N)		/*M�sN��̂Q�����z�������������֐�*/
//{
//
//	for (size_t i = 0; i < M; i++)
//	{
//		for (size_t j = 0; j < N; j++)
//		{
//			Mtr[i * N + j] = 0.0000;
//		}
//	}
//}
//
//void Make_H_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Hm)		/*H matrix�쐬�p�֐�*/
//{
//	//int i;	/*�s*/
//	//int j;	/*��*/
//
//	for (size_t i = 0; i < N_node; i++)
//	{
//		for (size_t j = 0; j < NUM_sub; j++)
//		{
//			//printf("\ni = %d\nj = %d\n",i,j);
//			if (i == j)
//			{
//				//printf("\ni = %d\nj = %d\n",i,j);
//				Hm[i + j * N_node] = -1.0;
//			}
//			else
//			{
//				Hm[i + j * N_node] = 0.0;
//			}
//		}
//	}
//	for (size_t i = 0; i < N_node; i++)
//	{
//		for (size_t j = NUM_sub; j < N_branch; j++)
//		{
//			if (i == branch[j].Node_pos)
//			{
//				Hm[i + j * N_node] = 1.0;
//			}
//			else if (i == branch[j].Node_neg)
//			{
//				Hm[i + j * N_node] = -1.0;
//			}
//			else
//			{
//				Hm[i + j * N_node] = 0.0;
//			}
//		}
//	}
//}
//
//void Make_Y_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Y)		/*Y matrix�쐬�p�֐�*/
//{
//	//int i;	/*�s*/
//	//int j;	/*��*/
//
//	for (size_t i = 0; i < N_branch; i++)
//	{
//		for (size_t j = 0; j < N_branch; j++)
//		{
//			//printf("\ni = %d\nj = %d\n",i,j);
//			if (i == j)
//			{
//				//printf("\ni = %d\nj = %d\n",i,j);
//				if (branch[j].r < 0.0001)
//				{
//					Y[i + j * N_branch] = 10000.0;
//				}
//				else
//				{
//					Y[i + j * N_branch] = 1.0 / branch[j].r;		/*�����ŃR���_�N�^���X[S]�ɕϊ�*/
//				}
//			}
//			else
//			{
//				Y[i + j * N_branch] = 0.0;
//			}
//		}
//	}
//}
//
//void Make_In_vector(const std::unique_ptr<NODE[]>& data, const std::unique_ptr<double[]>& In)		/*In vector�쐬�p�֐�*/
//{
//
//	for (size_t i = 0; i < N_node; i++)
//	{
//		In[i] = data[i].I;
//	}
//}
//
//void Transpose(const std::unique_ptr<double[]>& trans, const std::unique_ptr<double[]>& X, size_t row, size_t column)		/*�]�u�s��v�Z�p�֐�*/
//{
//	//int i, j;
//	double temp;
//
//	for (size_t i = 0; i < row; i++)
//	{
//		for (size_t j = 0; j < column; j++)
//		{
//			temp = X[i + j * row];
//			trans[j + i * column] = temp;
//		}
//	}
//}
//
////int mySimulate(std::vector<double> wait_time, FILE* file)
////{
////	/***  �e�s��錾��  ***/
////	double H[N_node * N_branch];		//Connection Matrix
////	double H_tra[N_branch * N_node];	//Transpose matrix of H
////	double Y[N_branch * N_branch];		//Conductance Matrix
////	double A[N_node * N_node];			//A matrix(H*Y*H_tra)
////	double A_inv[N_node * N_node];		//Inverse of A matrix
////
////	double TEMP[N_node * N_branch];	//Tast
////
////	double In[N_node];					//Node Current Matrix
////	double In_cpy[N_node];				//Copy of Node Current Matrix
////	double Vn[N_node];					//Node Voltage Matrix
////	double Ib[N_branch];				//Branch Current Matrix
////	double Vb[N_branch];				//Branch Voltage Matrix
////
////	int i, j, k, l;
////	k = 0;
////	l = 0;
////
////	int ERROR = 0;						//��H�v�Z�̍Čv�Z����p�t���O�iflag��S�ĉ��Z����NUM_sub�ȏ�Ȃ��H�v�Z���[�v�𔲂���BNUM_sub�ȉ��Ȃ�Čv�Z�B�j
////	int ERROR_REG = 0;
////	int count_loop = 0;
////	int count_loop_rec = 0;
////	//	int test_flag = 0;
////
////		/***  �}�g���b�N�X�v�Z�֌W  ***/
////	const char trans = 'N'; //Normal�ȏꍇ�B����͋K��l'N','T','C'����I��
////	int K_node = N_node;
////	int K_branch = N_branch;
////	int ld_node = N_node;
////	int ld_branch = N_branch;
////	double alpha = 1.0;
////	double beta = 0.0;
////
////	int* ipiv_Vn;
////	ipiv_Vn = (int*)calloc(N_node, sizeof(int));
////	int info_Vn;
////	int nrhs = 1;
////
////	int incx = 1;
////	int incy = 1;
////
////	/***  �\���̐錾��  ***/
////	TRAIN& tra[NUM_tra];				//�|�C���^�̔z��Ƃ��Ē�`
////	SUB& sub[NUM_sub];					//�|�C���^�̔z��Ƃ��Ē�`
////	NODE* node[N_node];
////	NODE* node_order[N_node];
////	BRANCH* branch[N_branch];
////	//	ESD *esd[NUM_esd];					//�|�C���^�̔z��Ƃ��Ē�`
////	//	INV *inv[NUM_inv];					//�|�C���^�̔z��Ƃ��Ē�`
////
////		/*** ���d��H���d�����i�[�ϐ��@***/
////	//double Lcir = 0.0;							//[J]
////	//double Lcir_c = 0.0;						//[kWh]
////
////
////	/*** ��ԍ\���̔z�񏉊����p�̍\���̔z��錾�� ***/
////	/*�I�t�s�[�N�_�C��*/
////	INI_TRA ini_tra[NUM_tra] = {	//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���, ��Ԏ��
////		//�R�����14�F45-15:30����ɂ���
////
////		//�����F�r�܁��V�h���a�J	
////		{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 0.0, 0.0, 1.0, 1080.0, 0, 0},				//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
////		{32.78 * 1000.0, 26.25 * 1000.0,  6.0, 4.0, 0.0, 0.0, 1.0, 1380.0, 0, 0},				//M�ԏd��[kg/��]�CT�ԏd��[kg/��]�CM�ԗ����CT�ԗ����C�����ʒu[m]�C�������x[km/h]�C�㉺����C���Ԏ���
////
////	};
////
////
////
////	/*** �ϓd���\���̔z�񏉊����p�̍\���̔z��錾 ***/
////	INI_SUB ini_sub[NUM_sub] = {
////		{-300.0, ess1, Iss_rated1, Vss1, 0.0, 1},						//[�V�S�����u]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d��, �ː��d���C�_�C�I�[�h
////		{5100.0, ess2, Iss_rated2, Vss2, 0.0, 1},						//[�i�R(��)]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d��, �ː��d���C�_�C�I�[�h
////		{8500.0, ess3, Iss_rated3, Vss3, 0.0, 1},					//[���ؓc]	�ϓd���ʒu�C�d���ϓ����C�����o�͒�R�C�o�͓d���C�o�͓d���C�ː��d���C�_�C�I�[�h
////	};
////
////	/***  �w�̍\���̔z��錾�� (�V�S�����u�w�) ***/
////	STATION sta[] = {	//����e�w
////		{"basement1",	-100.0},			//�w��, �w�ʒu	0
////		{"Ikebukuro",  	0.0},				//�w��, �w�ʒu	1
////		{"Mejiro",	    1200.0},	    	//�w��, �w�ʒu	2
////		{"Takadanobaba",  2100.0},			//�w��, �w�ʒu	3
////		{"Shin-Okubo",	3500.0},			//�w��, �w�ʒu	4
////		{"Shinjuku",	4800.0},			//�w��, �w�ʒu	5
////		{"Yoyogi",	5500.0},		    	//�w��, �w�ʒu	6
////		{"Harajuku",	7000.0},			//�w��, �w�ʒu	7
////		{"Shibuya",	8200.0},			    //�w��, �w�ʒu	8
////		{"basement2",	8300.0},			//�w��, �w�ʒu	9
////	};
////
////
////	sim_endflag = 0;
////
////	/***  CSV�t�@�C��1�I�[�v��&��s�ڏ�������  ***/
////	while (sim_endflag == 0) {
////
////		
////		/***  ��ԍ\���́E�ϓd���\���́EESD�\���̏�����  ***/
////		for (i = 0; i < NUM_sub; i++)
////		{
////			Make_substation(sub[i], &(ini_sub[i]), i);				//�i���o�����O�͒ʂ��ԍ�
////		}
////		for (i = 0; i < NUM_tra; i++)
////		{
////			Make_train(tra[i], &(ini_tra[i]), i + NUM_sub);							//�i���o�����O�͒ʂ��ԍ�
////		}
////
////		for (i = 0; i < NUM_tra; i++)
////		{
////			Calculate_BEmax(tra[i]);
////		}
////
////		////////////////////////// ���l�v�Z�� /////////////////////////
////		for (t = 0.0; t <= 60.0 * 45.0; t += dt)	/***** 1���ԃV�~�����[�V���� *****/
////		{
////			minute = t / 60.0;
////			count_loop_rec = 0;
////
////			////////////////////////// �ԗ��^���������v�Z�� /////////////////////
////				/*** �ԗ����s�p�^�[���͋[ & �^���������v�Z ***/
////				/*** �ԗ��p���[�v�Z ***/
////
////			for (i = 0; i < NUM_tra; i++)
////			{
////				Calculate_BEmax(tra[i]);
////
////				Run_pattern(tra[i], sta, t, wait_time);
////
////			}
////
////			/////////////////////////////////////////////////////////////////////
////
////		//	for (i = 0; i < NUM_sub; i++)
////		//	{
////		//		sub[i].Jss = sub[i].Jss_ini;
////		//	}
////
////			/*	if (t > 1020.05 && t < 1020.09) {
////					printf("-----------t = %f--------------\n", t);
////					printf("-----------accelflag vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%d\n", tra[i].accelflag);
////					}
////					printf("-----------Fmot vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%f\n", tra[i].Fmot);
////					}
////					printf("\n");
////					printf("\n");
////					printf("\n");
////				}*/
////
////			while (1)
////			{
////			start_while:				//�_�C�I�[�h���肪�Ԉ���Ă����ꍇ�C��������Čv�Z����
////
////				if (count_loop > 300000)
////				{
////					printf("\n");
////					printf("-----------t = %f--------------\n", t);
////					printf("-----------accelflag & position vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%d, %f\n", tra[i].accelflag, tra[i].x);
////					}
////					/*				printf("-----------Fmot vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%f\n", tra[i].Fmot);
////									}*/
////					printf("-----------Pmot vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%f\n", tra[i].Pmot);
////					}
////					printf("-----------Vp vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%f\n", tra[i].vp);
////					}
////					printf("-----------Vfc vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%f\n", tra[i].vfc);
////					}
////					printf("-----------ifl vehicle--------------\n");
////					for (i = 0; i < NUM_tra; i++)
////					{
////						printf("%f\n", tra[i].ifl);
////					}
////					/*		printf("-----------ifc vehicle--------------\n");
////							for (i = 0; i < NUM_tra; i++)
////							{
////								printf("%f\n", tra[i].ifc);
////							}
////							printf("-----------iinv vehicle--------------\n");
////							for (i = 0; i < NUM_tra; i++)
////							{
////								printf("%f\n", tra[i].iinv);
////							}*/
////							/*				printf("-----------vss substation--------------\n");
////											for (i = 0; i < NUM_sub; i++)
////											{
////												printf("%f\n", sub[i].vout);
////											}
////											printf("-----------Jss substation--------------\n");
////											for (i = 0; i < NUM_sub; i++)
////											{
////												printf("%f\n", sub[i].Jss);
////											}
////											printf("-----------Jss_ini substation--------------\n");
////											for (i = 0; i < NUM_sub; i++)
////											{
////												printf("%f\n", sub[i].Jss_ini);
////											}
////											printf("-----------iss substation--------------\n");
////											for (i = 0; i < NUM_sub; i++)
////											{
////												printf("%f\n", sub[i].iss);
////											}*/
////					printf("\n");
////					printf("\n");
////				}
////
////				if (count_loop > 300010) {
////					printf("Holy shit!!\n");
////					exit(EXIT_FAILURE);
////				}
////
////				/***  node[]�̏�����  ***/
////				for (i = 0; i < NUM_sub; i++)
////				{
////					Make_NODE_SS(node[i], sub[i]);
////					Make_NODE_SS(&(node_order[i]), sub[i]);
////				}
////
////				for (i = 0; i < NUM_tra; i++)
////				{
////					Make_NODE_TRAIN(&(node[i + NUM_sub]), tra[i]);
////					Make_NODE_TRAIN(&(node_order[i + NUM_sub]), tra[i]);
////				}
////
////				for (i = 0; i < N_branch; i++)
////				{
////					Initialize_BRANCH(&(branch[i]));
////				}
////
////				/*** �eMatrix�̏����� ***/
////				Initialize_matrix1(H, N_node, N_branch);
////				Initialize_matrix1(H_tra, N_branch, N_node);
////				Initialize_matrix1(Y, N_branch, N_branch);
////				Initialize_matrix1(A, N_node, N_node);
////				Initialize_matrix1(A_inv, N_node, N_node);
////				Initialize_matrix1(TEMP, N_node, N_branch);
////				Initialize_matrix1(In, N_node, 1);
////				Initialize_matrix1(In_cpy, N_node, 1);
////				Initialize_matrix1(Vn, N_node, 1);
////				Initialize_matrix1(Ib, N_branch, 1);
////				Initialize_matrix1(Vb, N_branch, 1);
////
////				/***  node_order[]���������Ƀ\�[�g  ***/
////				sort(node_order, N_node);
////
////				Make_branch(branch, node_order, tra);
////
////				//H�s��̐���(�u�񐶎Ԃ��܂ނ��d�V�X�e���̌���Ƃ�����vp.19�Q��)
////				Make_H_matrix(branch, H);
////
////				//�A�h�~�^���X�s��̐���
////				Make_Y_matrix(branch, Y);
////
////				//�m�[�h�d���s��̐���
////				Make_In_vector(node, In);
////				Make_In_vector(node, In_cpy);		//"dgetrs"�ł͉E�Ӄx�N�g�����������̉��ŏ㏑������Ă��܂��̂ŁC�v�Z�p��In���R�s�[����In_cpy���쐬
////
////				//H�s���]�u����
////				Transpose(H_tra, H, N_node, N_branch);
////
////				//H�s��~Y�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
////				dgemm(&trans, &trans, &K_node, &K_branch, &K_branch, &alpha, H, &ld_node, Y, &ld_branch, &beta, TEMP, &ld_node);
////
////				//��L�s��~�]�uH�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
////				dgemm(&trans, &trans, &K_node, &K_node, &K_branch, &alpha, TEMP, &ld_node, H_tra, &ld_branch, &beta, A, &ld_node);
////
////				//�A��1�����������������߂�H�s��~Y�s��~�]�uH�s��̉�A(N�sN��)��LU����
////				dgetrf(&ld_node, &K_branch, A, &ld_node, ipiv_Vn, &info_Vn);
////
////				//dgetrs�ŘA��1���������uA*Vn=In�v�������CIn�ɉ�(�m�[�h�d���s��)���㏑��
////				dgetrs(&trans, &K_node, &nrhs, A, &ld_node, ipiv_Vn, In_cpy, &ld_node, &info_Vn);	//In_cpy�ɕ������̉�Vn���i�[����Ă��邱�Ƃɒ���
////
////				Make_Y_matrix(branch, Y);		//�Ȃ���Vn�����߂�ۂ�"dgetrf"��Y-matrix�������������Ă��邽�߁C�����ōČv�Z�����Ă�i���ۂ͖��ʂȌv�Z�Ȃ̂ł������Ȃ����c�j
////
////				for (i = 0; i < N_node; i++)	/*Vn�x�N�g���̓}�C�i�X������*/
////				{
////					Vn[i] = -In_cpy[i];
////				}
////
////				//�]�uH�s��~Vn�x�N�g���Ńu�����`�d���x�N�g�����v�Z
////				dgemv(&trans, &K_branch, &K_node, &alpha, H_tra, &ld_branch, Vn, &incx, &beta, Vb, &incy);
////
////				//Y�s��~�u�����`�d���x�N�g���Ńu�����`�d���x�N�g�����v�Z
////				dgemv(&trans, &K_branch, &K_branch, &alpha, Y, &ld_branch, Vb, &incx, &beta, Ib, &incy);
////
////
////				for (i = 0; i < N_node; i++)			//�m�[�h�d���̌v�Z���ʂ��m�[�h�\���̂֕Ԃ��i���ʓI�ɂ́CVn�z��̓m�[�h�\���̂̏��Ԓʂ�ɕ���ł���H�j
////				{
////					node[i].V = Vn[i];
////				}
////
////				for (i = 0; i < NUM_sub; i++)			//�m�[�h�d���E�d���ƃu�����`�d������C�ϓd���o�͓d���E�d�����v�Z
////				{
////					sub[i].vout = node[i].V;
////					sub[i].iss = -In[i] + Ib[i];
////					Error_Detection_SS(sub[i]);
////					ERROR += sub[i].flag;
////				}
////
////				for (i = 0; i < NUM_tra; i++)
////				{
////					tra[i].vfc_old = tra[i].vfc;
////
////					for (j = 0; j < N_node; j++)
////					{
////						if (node[j]->Number == tra[i].name_train)	//�m�[�h�ԍ��Ɨ�Ԏ��ʔԍ�����v���Ă�����m�[�h�d�����ԃp���^�_�d���ɕԂ�
////						{
////							tra[i].vp = node[j]->V;
////							tra[i].vfc = tra[i].vp - tra[i].Rfl * tra[i].ifl;
////
////							//�����́E�u���[�L�́E���[�^�d���E���[�^�d�����v�Z
////							Calculation_traction_force_brake_force(tra[i]);
////
////							//dq���d���E�d�����烂�[�^�p���[���v�Z
////							Calculation_power_motor(tra[i]);
////
////							//���H�v�Z
////							Calculation_traction_circuit(tra[i]);
////						}
////					}
////				}
////
////				for (i = 0; i < NUM_tra; i++)
////				{
////					tra[i].judgement = fabs(tra[i].vfc - tra[i].vfc_old) / tra[i].vfc_old;
////
////					if (tra[i].judgement < 0.001) {
////						tra[i].reg_flag = 1;
////					}
////					else {
////						tra[i].reg_flag = 0;
////					}
////
////					if (tra[i].vfc > Vcmax) tra[i].reg_flag = 0;
////
////					ERROR_REG += tra[i].reg_flag;
////				}
////
////				/*		for (i = 0; i < NUM_tra; i++)
////						{
////							if (tra[i].accelflag == 3 || (tra[i].accelflag == 5 && tra[i].v_c > tra[i].con_speed)) {
////								Error_Detection_REG(tra[i]);
////							}
////							else {
////								tra[i].reg_flag = 1;
////							}
////
////							ERROR_REG += tra[i].reg_flag;
////						}*/
////
////						/*		if (t > 1020.05 && t < 1020.09) {
////									printf("-----------t = %f--------------\n", t);
////									printf("-----------accelflag vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%d\n", tra[i].accelflag);
////									}
////									printf("-----------Fmot vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%f\n", tra[i].Fmot);
////									}
////									printf("-----------Pmot vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%f\n", tra[i].Pmot);
////									}
////									printf("-----------Vp vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%f\n", tra[i].vp);
////									}
////									printf("-----------ifl vehicle--------------\n");
////									for (i = 0; i < NUM_tra; i++)
////									{
////										printf("%f\n", tra[i].ifl);
////									}
////									printf("-----------vss substation--------------\n");
////									for (i = 0; i < NUM_sub; i++)
////									{
////										printf("%f\n", sub[i].vout);
////									}
////									printf("-----------iss substation--------------\n");
////									for (i = 0; i < NUM_sub; i++)
////									{
////										printf("%f\n", sub[i].iss);
////									}
////									printf("\n");
////									printf("\n");
////									printf("\n");
////								}*/
////
////				if (ERROR < NUM_sub || ERROR_REG < NUM_tra)
////				{
////					for (i = 0; i < N_node; i++)
////					{
////						free(node[i]);
////						free(node_order[i]);
////					}
////					for (i = 0; i < N_branch; i++)
////					{
////						free(branch[i]);
////					}
////
////					ERROR = 0;
////					ERROR_REG = 0;
////					count_loop++;
////					goto start_while;					//while���̐擪�ɖ߂�
////				}
////				else
////				{
////					ERROR = 0;
////					ERROR_REG = 0;
////					count_loop_rec = count_loop;
////					count_loop = 0;
////					break;
////				}
////			}
////
////			/***  ���̏�Ԃ��v�Z  ***/
////			for (i = 0; i < NUM_tra; i++)
////			{
////				if (tra[i].accelflag != 4) Solve_next_state(tra[i]);
////			}
////
////			/***  �ԗ��G�l���M�[�v�Z  ***/
////			for (i = 0; i < NUM_tra; i++)
////			{
////				Calculation_power_train(tra[i]);
////				Calculation_energy_train(tra[i]);
////			}
////
////			for (i = 0; i < NUM_sub; i++)
////			{
////				Calculation_power_sub(sub[i]);
////				Calculation_energy_sub(sub[i]);
////			}
////
////			for (i = 0; i < N_node; i++)
////			{
////				free(node[i]);
////				free(node_order[i]);
////			}
////
////			for (i = 0; i < N_branch; i++)
////			{
////				free(branch[i]);
////			}
////
////
////			////////////////////////// CSV�t�@�C���o�͕� /////////////////////////
////
////
////			//////////////////////////////////////////////////////////////////////
////
////
////			/////////////////////////////  cmd�\����  ////////////////////////////
////
////			//////////////////////////////////////////////////////////////////////
////		}
////		/********** ���l�v�Z���[�v�I�� **********/
////
////
////		double sub_total = 0.0;
////
////		for (i = 0; i < NUM_sub; i++)
////		{
////			sub_total += sub[i].Ess_c;
////		}
////
////		if (Consider_loop_flag == 0) {
////			fprintf(file, "%.1f, %f, %f, %f, %f, %f\n", wait_time[0], wait_time[1], wait_time[2], wait_time[3], wait_time[4], sub_total);
////		}
////
////
////		/*****  �������J��  *****/
////		for (i = 0; i < NUM_tra; i++)
////		{
////			free(tra[i]);
////		}
////		for (i = 0; i < NUM_sub; i++)
////		{
////			free(sub[i]);
////		}
////
////
////
////		sim_endflag = 1;
////
////		
////	}
////
////
////	return  0;
////
////}
//
//constexpr char trans = 'N'; //Normal�ȏꍇ�B����͋K��l'N','T','C'����I��
//
//int mySimulate(std::vector<double> wait_time, CsvWriter& _cw, std::mutex& _mtx)
//{
//	/***  �e�s��錾��  ***/
//	/***  �e�s��錾��  ***/
//	std::unique_ptr<double[]> H = std::make_unique<double[]>(N_node * N_branch);		//Connection Matrix
//	std::unique_ptr<double[]> H_tra = std::make_unique<double[]>(N_branch * N_node);	//Transpose matrix of H
//	std::unique_ptr<double[]> Y = std::make_unique<double[]>(N_branch * N_branch);		//Conductance Matrix
//	std::unique_ptr<double[]> A = std::make_unique<double[]>(N_node * N_node);			//A matrix(H*Y*H_tra)
//	//std::unique_ptr<double[]> A_inv[N_node * N_node);		//Inverse of A matrix
//
//	std::unique_ptr<double[]> TEMP = std::make_unique<double[]>(N_node * N_branch);	//Test
//
//	std::unique_ptr<double[]> In = std::make_unique<double[]>(N_node);					//Node Current Matrix
//	std::unique_ptr<double[]> In_cpy = std::make_unique<double[]>(N_node);				//Copy of Node Current Matrix
//	std::unique_ptr<double[]> Vn = std::make_unique<double[]>(N_node);					//Node Voltage Matrix
//	std::unique_ptr<double[]> Ib = std::make_unique<double[]>(N_branch);				//Branch Current Matrix
//	std::unique_ptr<double[]> Vb = std::make_unique<double[]>(N_branch);
//
//	int i = 0, j = 0;
//
//	int ERROR = 0;						//��H�v�Z�̍Čv�Z����p�t���O�iflag��S�ĉ��Z����NUM_sub�ȏ�Ȃ��H�v�Z���[�v�𔲂���BNUM_sub�ȉ��Ȃ�Čv�Z�B�j
//	int ERROR_REG = 0;
//	int count_loop = 0;
//	int count_loop_rec = 0;
//	//	int test_flag = 0;
//
//		/***  �}�g���b�N�X�v�Z�֌W  ***/
//	
//	int K_node = N_node;
//	int K_branch = N_branch;
//	int ld_node = N_node;
//	int ld_branch = N_branch;
//	double alpha = 1.0;
//	double beta = 0.0;
//
//	//int* ipiv_Vn;
//	std::unique_ptr<int[]> ipiv_Vn = std::make_unique<int[]>(N_node);
//	int info_Vn;
//	int nrhs = 1;
//
//	int incx = 1;
//	int incy = 1;
//
//	/***  �\���̐錾��  ***/
//	std::unique_ptr<TRAIN[]> tra = std::make_unique<TRAIN[]>(NUM_tra);	//�|�C���^�̔z��Ƃ��Ē�`
//	std::unique_ptr<SUB[]> sub = std::make_unique<SUB[]>(NUM_sub);
//	std::unique_ptr<NODE[]> node = std::make_unique<NODE[]>(N_node);
//	std::unique_ptr<NODE[]> node_order = std::make_unique<NODE[]>(N_node);
//	std::unique_ptr<BRANCH[]> branch = std::make_unique<BRANCH[]>(N_branch);
//	std::unique_ptr<ESD[]> esd = std::make_unique<ESD[]>(NUM_esd);
//	//	ESD *esd[NUM_esd];					//�|�C���^�̔z��Ƃ��Ē�`
//	//	INV *inv[NUM_inv];					//�|�C���^�̔z��Ƃ��Ē�`
//
//		/*** ���d��H���d�����i�[�ϐ��@***/
//	//double Lcir = 0.0;							//[J]
//	//double Lcir_c = 0.0;						//[kWh]
//
//
//	
//
//	sim_endflag = 0;
//
//	/***  CSV�t�@�C��1�I�[�v��&��s�ڏ�������  ***/
//	while (sim_endflag == 0) {
//
//
//		/***  ��ԍ\���́E�ϓd���\���́EESD�\���̏�����  ***/
//		for (i = 0; i < NUM_sub; i++)
//		{
//			Make_substation(sub[i], &(ini_sub[i]), i);				//�i���o�����O�͒ʂ��ԍ�
//		}
//		for (i = 0; i < NUM_tra; i++)
//		{
//			Make_train(tra[i], &(ini_tra[i]), i + NUM_sub);							//�i���o�����O�͒ʂ��ԍ�
//		}
//
//		for (i = 0; i < NUM_tra; i++)
//		{
//			Calculate_BEmax(tra[i]);
//		}
//
//		////////////////////////// ���l�v�Z�� /////////////////////////
//		for (t = 0.0; t <= 60.0 * 45.0; t += dt)	/***** 1���ԃV�~�����[�V���� *****/
//		{
//			minute = t / 60.0;
//			count_loop_rec = 0;
//
//			////////////////////////// �ԗ��^���������v�Z�� /////////////////////
//				/*** �ԗ����s�p�^�[���͋[ & �^���������v�Z ***/
//				/*** �ԗ��p���[�v�Z ***/
//
//			for (i = 0; i < NUM_tra; i++)
//			{
//				Calculate_BEmax(tra[i]);
//
//				Run_pattern(tra[i], sta, t, wait_time);
//
//			}
//
//			/////////////////////////////////////////////////////////////////////
//
//		//	for (i = 0; i < NUM_sub; i++)
//		//	{
//		//		sub[i].Jss = sub[i].Jss_ini;
//		//	}
//
//			/*	if (t > 1020.05 && t < 1020.09) {
//					printf("-----------t = %f--------------\n", t);
//					printf("-----------accelflag vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%d\n", tra[i].accelflag);
//					}
//					printf("-----------Fmot vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%f\n", tra[i].Fmot);
//					}
//					printf("\n");
//					printf("\n");
//					printf("\n");
//				}*/
//
//			while (1)
//			{
//			start_while:				//�_�C�I�[�h���肪�Ԉ���Ă����ꍇ�C��������Čv�Z����
//
//				if (count_loop > 300000)
//				{
//					printf("\n");
//					printf("-----------t = %f--------------\n", t);
//					printf("-----------accelflag & position vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%d, %f\n", tra[i].accelflag, tra[i].x);
//					}
//					/*				printf("-----------Fmot vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%f\n", tra[i].Fmot);
//									}*/
//					printf("-----------Pmot vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%f\n", tra[i].Pmot);
//					}
//					printf("-----------Vp vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%f\n", tra[i].vp);
//					}
//					printf("-----------Vfc vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%f\n", tra[i].vfc);
//					}
//					printf("-----------ifl vehicle--------------\n");
//					for (i = 0; i < NUM_tra; i++)
//					{
//						printf("%f\n", tra[i].ifl);
//					}
//					/*		printf("-----------ifc vehicle--------------\n");
//							for (i = 0; i < NUM_tra; i++)
//							{
//								printf("%f\n", tra[i].ifc);
//							}
//							printf("-----------iinv vehicle--------------\n");
//							for (i = 0; i < NUM_tra; i++)
//							{
//								printf("%f\n", tra[i].iinv);
//							}*/
//							/*				printf("-----------vss substation--------------\n");
//											for (i = 0; i < NUM_sub; i++)
//											{
//												printf("%f\n", sub[i].vout);
//											}
//											printf("-----------Jss substation--------------\n");
//											for (i = 0; i < NUM_sub; i++)
//											{
//												printf("%f\n", sub[i].Jss);
//											}
//											printf("-----------Jss_ini substation--------------\n");
//											for (i = 0; i < NUM_sub; i++)
//											{
//												printf("%f\n", sub[i].Jss_ini);
//											}
//											printf("-----------iss substation--------------\n");
//											for (i = 0; i < NUM_sub; i++)
//											{
//												printf("%f\n", sub[i].iss);
//											}*/
//					printf("\n");
//					printf("\n");
//				}
//
//				if (count_loop > 300010) {
//					printf("Holy shit!!\n");
//					exit(EXIT_FAILURE);
//				}
//
//				/***  node[]�̏�����  ***/
//				for (i = 0; i < NUM_sub; i++)
//				{
//					Make_NODE_SS(node[i], sub[i]);
//					Make_NODE_SS(node_order[i], sub[i]);
//				}
//
//				for (i = 0; i < NUM_tra; i++)
//				{
//					Make_NODE_TRAIN(node[i + NUM_sub], tra[i]);
//					Make_NODE_TRAIN(node_order[i + NUM_sub], tra[i]);
//				}
//
//				for (i = 0; i < N_branch; i++)
//				{
//					Initialize_BRANCH(branch[i]);
//				}
//
//				/*** �eMatrix�̏����� ***/
//				Initialize_matrix1(H, N_node, N_branch);
//				Initialize_matrix1(H_tra, N_branch, N_node);
//				Initialize_matrix1(Y, N_branch, N_branch);
//				Initialize_matrix1(A, N_node, N_node);
//				Initialize_matrix1(TEMP, N_node, N_branch);
//				Initialize_matrix1(In, N_node, 1);
//				Initialize_matrix1(In_cpy, N_node, 1);
//				Initialize_matrix1(Vn, N_node, 1);
//				Initialize_matrix1(Ib, N_branch, 1);
//				Initialize_matrix1(Vb, N_branch, 1);
//
//				/***  node_order[]���������Ƀ\�[�g  ***/
//				//sort_s(node_order, N_node);
//				sort(node_order, N_node);
//
//				Make_branch(branch, node_order, tra);
//
//				//H�s��̐���(�u�񐶎Ԃ��܂ނ��d�V�X�e���̌���Ƃ�����vp.19�Q��)
//				Make_H_matrix(branch, H);
//
//				//�A�h�~�^���X�s��̐���
//				Make_Y_matrix(branch, Y);
//
//				//�m�[�h�d���s��̐���
//				Make_In_vector(node, In);
//				Make_In_vector(node, In_cpy);		//"dgetrs"�ł͉E�Ӄx�N�g�����������̉��ŏ㏑������Ă��܂��̂ŁC�v�Z�p��In���R�s�[����In_cpy���쐬
//
//				//H�s���]�u����
//				Transpose(H_tra, H, N_node, N_branch);
//
//				//H�s��~Y�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
//				dgemm(&trans, &trans, &K_node, &K_branch, &K_branch, &alpha, H.get(), &ld_node, Y.get(), &ld_branch, &beta, TEMP.get(), &ld_node);
//
//				//��L�s��~�]�uH�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
//				dgemm(&trans, &trans, &K_node, &K_node, &K_branch, &alpha, TEMP.get(), &ld_node, H_tra.get(), &ld_branch, &beta, A.get(), &ld_node);
//
//				//�A��1�����������������߂�H�s��~Y�s��~�]�uH�s��̉�A(N�sN��)��LU����
//				dgetrf(&ld_node, &K_branch, A.get(), &ld_node, ipiv_Vn.get(), &info_Vn);
//
//				//dgetrs�ŘA��1���������uA*Vn=In�v�������CIn�ɉ�(�m�[�h�d���s��)���㏑��
//				dgetrs(&trans, &K_node, &nrhs, A.get(), &ld_node, ipiv_Vn.get(), In_cpy.get(), &ld_node, &info_Vn);	//In_cpy�ɕ������̉�Vn���i�[����Ă��邱�Ƃɒ���
//
//				Make_Y_matrix(branch, Y);		//�Ȃ���Vn�����߂�ۂ�"dgetrf"��Y-matrix�������������Ă��邽�߁C�����ōČv�Z�����Ă�i���ۂ͖��ʂȌv�Z�Ȃ̂ł������Ȃ����c�j
//
//				for (i = 0; i < N_node; i++)	/*Vn�x�N�g���̓}�C�i�X������*/
//				{
//					Vn[i] = -In_cpy[i];
//				}
//
//				//�]�uH�s��~Vn�x�N�g���Ńu�����`�d���x�N�g�����v�Z
//				dgemv(&trans, &K_branch, &K_node, &alpha, H_tra.get(), &ld_branch, Vn.get(), &incx, &beta, Vb.get(), &incy);
//
//				//Y�s��~�u�����`�d���x�N�g���Ńu�����`�d���x�N�g�����v�Z
//				dgemv(&trans, &K_branch, &K_branch, &alpha, Y.get(), &ld_branch, Vb.get(), &incx, &beta, Ib.get(), &incy);
//
//
//				for (i = 0; i < N_node; i++)			//�m�[�h�d���̌v�Z���ʂ��m�[�h�\���̂֕Ԃ��i���ʓI�ɂ́CVn�z��̓m�[�h�\���̂̏��Ԓʂ�ɕ���ł���H�j
//				{
//					node[i].V = Vn[i];
//				}
//
//				for (i = 0; i < NUM_sub; i++)			//�m�[�h�d���E�d���ƃu�����`�d������C�ϓd���o�͓d���E�d�����v�Z
//				{
//					sub[i].vout = node[i].V;
//					sub[i].iss = -In[i] + Ib[i];
//					Error_Detection_SS(sub[i]);
//					ERROR += sub[i].flag;
//				}
//
//				for (i = 0; i < NUM_tra; i++)
//				{
//					tra[i].vfc_old = tra[i].vfc;
//
//					for (j = 0; j < N_node; j++)
//					{
//						if (node[j].Number == tra[i].name_train)	//�m�[�h�ԍ��Ɨ�Ԏ��ʔԍ�����v���Ă�����m�[�h�d�����ԃp���^�_�d���ɕԂ�
//						{
//							tra[i].vp = node[j].V;
//							tra[i].vfc = tra[i].vp - tra[i].Rfl * tra[i].ifl;
//
//							//�����́E�u���[�L�́E���[�^�d���E���[�^�d�����v�Z
//							Calculation_traction_force_brake_force(tra[i]);
//
//							//dq���d���E�d�����烂�[�^�p���[���v�Z
//							Calculation_power_motor(tra[i]);
//
//							//���H�v�Z
//							Calculation_traction_circuit(tra[i]);
//						}
//					}
//				}
//
//				for (i = 0; i < NUM_tra; i++)
//				{
//					tra[i].judgement = fabs(tra[i].vfc - tra[i].vfc_old) / tra[i].vfc_old;
//
//					if (tra[i].judgement < 0.001) {
//						tra[i].reg_flag = 1;
//					}
//					else {
//						tra[i].reg_flag = 0;
//					}
//
//					if (tra[i].vfc > Vcmax) tra[i].reg_flag = 0;
//
//					ERROR_REG += tra[i].reg_flag;
//				}
//
//				/*		for (i = 0; i < NUM_tra; i++)
//						{
//							if (tra[i].accelflag == 3 || (tra[i].accelflag == 5 && tra[i].v_c > tra[i].con_speed)) {
//								Error_Detection_REG(tra[i]);
//							}
//							else {
//								tra[i].reg_flag = 1;
//							}
//
//							ERROR_REG += tra[i].reg_flag;
//						}*/
//
//						/*		if (t > 1020.05 && t < 1020.09) {
//									printf("-----------t = %f--------------\n", t);
//									printf("-----------accelflag vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%d\n", tra[i].accelflag);
//									}
//									printf("-----------Fmot vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%f\n", tra[i].Fmot);
//									}
//									printf("-----------Pmot vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%f\n", tra[i].Pmot);
//									}
//									printf("-----------Vp vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%f\n", tra[i].vp);
//									}
//									printf("-----------ifl vehicle--------------\n");
//									for (i = 0; i < NUM_tra; i++)
//									{
//										printf("%f\n", tra[i].ifl);
//									}
//									printf("-----------vss substation--------------\n");
//									for (i = 0; i < NUM_sub; i++)
//									{
//										printf("%f\n", sub[i].vout);
//									}
//									printf("-----------iss substation--------------\n");
//									for (i = 0; i < NUM_sub; i++)
//									{
//										printf("%f\n", sub[i].iss);
//									}
//									printf("\n");
//									printf("\n");
//									printf("\n");
//								}*/
//
//				if (ERROR < NUM_sub || ERROR_REG < NUM_tra)
//				{
//					/*for (i = 0; i < N_node; i++)
//					{
//						free(node[i]);
//						free(node_order[i]);
//					}
//					for (i = 0; i < N_branch; i++)
//					{
//						free(branch[i]);
//					}*/
//
//					/*node = std::make_unique<NODE[]>(N_node);
//					node_order = std::make_unique<NODE[]>(N_node);
//					branch = std::make_unique<BRANCH[]>(N_branch);*/
//
//					ERROR = 0;
//					ERROR_REG = 0;
//					count_loop++;
//					//goto start_while;					//while���̐擪�ɖ߂�
//					continue;
//				}
//				else
//				{
//					ERROR = 0;
//					ERROR_REG = 0;
//					count_loop_rec = count_loop;
//					count_loop = 0;
//					break;
//				}
//			}
//
//			/***  ���̏�Ԃ��v�Z  ***/
//			for (i = 0; i < NUM_tra; i++)
//			{
//				if (tra[i].accelflag != 4) Solve_next_state(tra[i]);
//			}
//
//			/***  �ԗ��G�l���M�[�v�Z  ***/
//			for (i = 0; i < NUM_tra; i++)
//			{
//				Calculation_power_train(tra[i]);
//				Calculation_energy_train(tra[i]);
//			}
//
//			for (i = 0; i < NUM_sub; i++)
//			{
//				Calculation_power_sub(sub[i]);
//				Calculation_energy_sub(sub[i]);
//			}
//
//			/*for (i = 0; i < N_node; i++)
//			{
//				free(node[i]);
//				free(node_order[i]);
//			}
//
//			for (i = 0; i < N_branch; i++)
//			{
//				free(branch[i]);
//			}*/
//
//			/*node = std::make_unique<NODE[]>(N_node);
//			node_order = std::make_unique<NODE[]>(N_node);
//			branch = std::make_unique<BRANCH[]>(N_branch);*/
//
//			////////////////////////// CSV�t�@�C���o�͕� /////////////////////////
//
//
//			//////////////////////////////////////////////////////////////////////
//
//
//			/////////////////////////////  cmd�\����  ////////////////////////////
//
//			//////////////////////////////////////////////////////////////////////
//		}
//		/********** ���l�v�Z���[�v�I�� **********/
//
//
//		double sub_total = 0.0;
//
//		for (i = 0; i < NUM_sub; i++)
//		{
//			sub_total += sub[i].Ess_c;
//		}
//
//		if (Consider_loop_flag == 0) {
//			#define tos std::to_string
//			//fprintf(file, "%.1f, %f, %f, %f, %f, %f\n", wait_time[0], wait_time[1], wait_time[2], wait_time[3], wait_time[4], sub_total);
//			std::string wl(tos(wait_time[0]) + "," + tos(wait_time[1]) + "," + tos(wait_time[2]) + "," + tos(wait_time[3]) + "," + tos(wait_time[4]) + "," + tos(sub_total));
//			//std::lock_guard<std::mutex> _lock(_mtx);
//			_cw.PushTask(wl, _mtx);
//		}
//
//
//		/*****  �������J��  *****/
//		/*for (i = 0; i < NUM_tra; i++)
//		{
//			free(tra[i]);
//		}
//		for (i = 0; i < NUM_sub; i++)
//		{
//			free(sub[i]);
//		}*/
//
//		//std::unique_ptr<TRAIN[]> tra = std::make_unique<TRAIN[]>(NUM_tra);	//�|�C���^�̔z��Ƃ��Ē�`
//		//std::unique_ptr<SUB[]> sub = std::make_unique<SUB[]>(NUM_sub);
//
//		sim_endflag = 1;
//
//
//	}
//
//
//	return  0;
//
//}