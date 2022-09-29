#pragma once
#include "particle.h"

class scinti
{
    double limtozero(double num);
    std::vector<double> scintibc();
    bool initpointcheck(particle initptcl, particle ptcl);
    public:
        double pt_x_,
            pt_y_,
            pt_z_,
            dir_theta_,
            dir_phi_,
            depth_,
            rad_,
            z_,
            dens_,
            ndens_,
            atomweight_;
        std::vector<std::vector<double>> crosssec_table_;
        void initcs(std::string conffilepath);
        double crosssec(double ene, int type);
        void initscinti(double pt_x, double pt_y, double pt_z, double theta, double phi, double depth, double z, double dens, double atomweight);
        std::vector<std::vector<double>> intersec(particle ptcl);
        double intersec_dist(particle initptcl, particle ptcl);
};