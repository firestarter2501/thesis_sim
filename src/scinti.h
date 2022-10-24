#pragma once
#include "particle.h"
#include "simfunc.h"

class scinti
{
    std::vector<std::vector<double>> intersec(particle ptcl);
    std::vector<bool> intersecenablecheck(std::vector<std::vector<double>> intersecvec, particle ptcl);
    std::string showfacetype(int type);
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
            atomweight_,
            pmtsdevslope,
            pmtsdevintersec,
            ene_buffer_;
        std::vector<std::vector<double>> crosssec_table_;
        void initscinti(std::string data, std::string csfile);
        void addscintidata(double pt_x, double pt_y, double pt_z, double theta, double phi, double depth, double z, double dens, double atomweight, double pmtsdevslope, double pmtsdevintersec);
        double intersec_dist(particle ptcl);
        double ptclinsidecheck(particle ptcl);
        double scintillation(particle &ptcl);
};