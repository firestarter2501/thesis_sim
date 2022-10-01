#pragma once
#include <iostream>
#include <fstream>
#include <random>
#include <Eigen/Dense>

class particle
{
    double limtozero(double num);
    public:
        double ene_,
            pt_x_,
            pt_y_,
            pt_z_,
            dir_theta_,
            dir_phi_;
        void initptcl(double ene, double x, double y, double z);
        void initptcl_test(int num);
        void move(double dist);
        void move_test(double dist);
        void turn(double angle);
        void turn_test(double theta, double phi, double angle);
};