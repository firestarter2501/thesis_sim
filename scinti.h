#pragma once
#include "particle.h"
#include <list>

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
        list<double> intersec(particle ptcl);
};

