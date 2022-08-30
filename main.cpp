#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <omp.h>
#include "simfunc.h"
#include "particle.h"
#include "scinti.h"

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

int main()
{
    std::vector<particle> ray_list;
    std::vector<scinti> scintillator;
    std::vector<double> initptclconf_list;
    std::vector<double> initscinticonf_list;
    std::ifstream initptclconf("initptcl.conf");
    std::ifstream initscinticonf("initscinti.conf");
    std::string initline;
    while (std::getline(initptclconf, initline))
    {
        initptclconf_list.push_back(std::stod(initline));
        std::cout << initline << std::endl;
    }
    for(int num = 0; num < initptclconf_list.size()/4; num++)
    {
        ray_list.push_back(particle());
        ray_list.back().initptcl(initptclconf_list.at(0+(num*4)), initptclconf_list.at(1+(num*4)), initptclconf_list.at(2+(num*4)), initptclconf_list.at(3+(num*4)));
    }
    std::cout << "---" << std::endl;
    while (std::getline(initscinticonf, initline))
    {
        initscinticonf_list.push_back(std::stod(initline));
        std::cout << initline << std::endl;
    }
    for(int num = 0; num < initscinticonf_list.size()/9; num++)
    {
        scintillator.push_back(scinti());
        scintillator.back().initscinti(initscinticonf_list.at(0+(num*9)), initscinticonf_list.at(1+(num*9)), initscinticonf_list.at(2+(num*9)), initscinticonf_list.at(3+(num*9))*M_PI, initscinticonf_list.at(4+(num*9))*M_PI, initscinticonf_list.at(5+(num*9)), initscinticonf_list.at(6+(num*9)), initscinticonf_list.at(7+(num*9)), initscinticonf_list.at(8+(num*9)));
    }
    std::cout << "---" << std::endl;
    std::cout << scintillator.size() << std::endl;

    // for (int run = 0; run < run_count; run++)
    // {
    //     std::vector<particle> photon;
    //     photon.push_back(ray_list.at(0));
    //     photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
    //     while (0 < photon.back().ene_)
    //     {
    //         double total_traject_dist = 0;
    //         for (int scinti_num = 0; scinti_num < scintillator.size(); scinti_num++)
    //         {
    //             double traject_dist = scintillator.at(scinti_num).intersec_dist(photon.back());
    //             total_traject_dist += traject_dist;
    //             if (traject_dist < 0.01)
    //             {
    //                 continue;
    //             }
    //             double pe_cs = pe_crosssec(photon.back().ene_, scintillator.at(scinti_num).z_),
    //                 cs_ang = cs_angle(photon.back().ene_),
    //                 cs_cs = kleinnishinaeq(photon.back().ene_, cs_ang);
    //         }
    //         if (total_traject_dist < 0.01)
    //         {
    //             break;
    //         }
    //     }
    // }
}