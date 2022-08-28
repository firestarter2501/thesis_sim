#pragma once

// 光電吸収断面積を返す関数
double pe_crosssec(double ene, int z);

// コンプトン散乱するときの散乱光子のエネルギーを返す関数
double scatphotonene(double ene, double angle);

// クライン仁科の式を計算する関数
double kleinnishinaeq(double ene, double angle);

// von Neumannの棄却法でコンプトン散乱の角度分布に基づいた散乱角の乱数を返す関数
double cs_angle(double ene);

void cs_angle_test(int num);

double reactlen(double crosssec, double dens);