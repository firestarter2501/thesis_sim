#pragma once
#include "particle.h"

class scinti
{
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
        void initsinti(double x, double y, double z, double theta, double phi, double depth, double dens, double atomweight);
        std::vector<std::vector<double>> intersec(particle ptcl);
        double intersec_dist(particle ptcl);
        bool internal_judge(particle ptcl);
};

