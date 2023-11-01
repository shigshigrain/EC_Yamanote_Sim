
#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <thread>



const std::vector<std::vector<double>> shift_wait = {

	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w1�ł̒������Ԉꗗ 
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w2�ł̒������Ԉꗗ
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w3
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w4
	{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w5
	//{-5.0, -2.5, 0.0, 2.5 , 5.0}, // �w6�ȍ~���L�q�\�@�c�����̐��͉w�̐��ƈ�v

};


static std::queue<std::vector<double>> wait_que;

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

	// �V�~�����[�V�����̒��g
	// �����ŁA�w�ł̒������ԃ^�C�~���O�͈ȉ��Ŏ擾�ł���

	// main�֐��ɏ����Ă��������������ӂɏ�肭�\��t����B������Ȃ������瑊�k

	std::vector<double> wait_time = wait_que.front(); // �擪�̃f�[�^��ǂݎ��B���̃��[�v�ł͐擪�̃f�[�^���̂Ă邱�ƂŎ��̃f�[�^��ǂ݂�����
	// �����œǂ݂������e�w���Ƃ̒�����������肭�V�~�����[�V�����ɓK��������


	// ���ʂ��t�H���_�ɏo�͂���  �t�@�C���������j�[�N�ɕt���Ă���

	// �V�~�����[�V�����{�̂�main�֐��ɋL�q�����֐��ɂ��Ă����ƁA��������񉻂��邱�Ƃ��ł���悤�ɂȂ�(�����ł͂��Ă��Ȃ�)

	return;

}


// main�֐����ǂ�
void _main() {

	make_que_wait(); // ���������̗񋓃��X�g����

	std::thread pallarel;

	while (!wait_que.empty()) // ���X�g�ɃV�~�����[�V�������ׂ��������炵�p�^�[�����c���Ă���΁A���[�v���s
	{
		simulate(wait_que.front()); // ����wait_que�̐擪�f�[�^���g�p���Ĉ�񕪂̃V�~�����[�V�������I��
		
		//std::vector<double> temp = wait_que.pop(); // ���V�~�����[�V�������I������������炵�p�^�[���͔j���B���̃p�^�[��������ɐ擪�ɂȂ�



	}

	return;

}

