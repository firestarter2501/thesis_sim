#include "scinti.h"

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

using namespace std;
using namespace Eigen;

void scinti::initsinti(double x, double y, double z, double theta, double phi, double depth, double dens, double atomweight)
{
    this->pt_x_ = x;
    this->pt_y_ = y;
    this->pt_z_ = z;
    this->dir_theta_ = theta;
    this->dir_phi_ = phi;
    this->depth_ = depth;
    this->rad_ = depth / 2;
    this->z_ = z;
    this->dens_ = dens;
    this->ndens_ = (dens / atomweight) * NUMA;
    this->atomweight_ = atomweight;
}

//std::vector<double> scinti::intersec(particle ptcl)
//{
//    // 直線の定義
//    // 手前の判定面の定義と判定
//    // 奥の判定面の定義と判定
//    // 円筒面の定義と判定
//}