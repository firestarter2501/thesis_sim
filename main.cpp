#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <math.h>
#include "particle.h"

using namespace std;
using namespace Eigen;

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

int main()
{
    particle test = { 0, 0, 0, 0, 0, 0 };
    test.move_test(5);
}