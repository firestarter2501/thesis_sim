#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <omp.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "simfunc.h"
#include "particle.h"
#include "scinti.h"

using namespace std;
using namespace Eigen;

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

int main()
{
    int run_count = 10;
    vector<particle> ray_list;
    vector<scinti> scintillator;

    ray_list.push_back(particle());
    ray_list.back().initptcl(661.6, 5, 0, 0);

    for (int run = 0; run < run_count; run++)
    {
        vector<particle> photon;
        photon.push_back(ray_list.at(0));
        photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
        while (0 < photon.back().ene_)
        {
            for (int scinti_num = 0; scinti_num < scintillator.size(); scinti_num++)
            {

            }
        }
    }
}