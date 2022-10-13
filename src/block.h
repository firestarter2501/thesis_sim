#pragma once
#include "particle.h"
#include "simfunc.h"

class block
{
    public:
        double pt_x_,
            pt_y_,
            pt_z_,
            height_,
            width_,
            depth_,
            z_,
            dens_,
            ndens_,
            atomweight_;
        std::vector<std::vector<double>> crosssec_table_;
        void initcs(std::string conffilepath);
        double crosssec(double ene, int type);
        void initblock(double pt_x, double pt_y, double pt_z, double height, double width, double depth, double z, double dens, double atomweight);
        double intersec_dist(particle ptcl);
        double ptclinsidecheck(particle ptcl);
        void react(std::string initblock, std::string initcs, particle &ptcl, bool &react_flag, bool &absorp_flag);
};