#pragma once

// ���d�z���f�ʐς�Ԃ��֐�
double pe_crosssec(double ene, int z);

// �R���v�g���U������Ƃ��̎U�����q�̃G�l���M�[��Ԃ��֐�
double scatphotonene(double ene, double angle);

// �N���C���m�Ȃ̎����v�Z����֐�
double kleinnishinaeq(double ene, double angle);

// von Neumann�̊��p�@�ŃR���v�g���U���̊p�x���z�Ɋ�Â����U���p�̗�����Ԃ��֐�
double cs_angle(double ene);

void cs_angle_test(int num);

double reactlen(double crosssec, double dens);