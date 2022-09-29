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

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

int main()
{
    std::vector<particle> ray_list;
    std::vector<scinti> scintillator;
    std::vector<double> initptclconf_list;
    std::vector<double> initscinticonf_list;
    std::ifstream initptclconf("./data/initptcl.conf");
    std::ifstream initscinticonf("./data/initscinti.conf");
    std::string initline;
    int run_count = 0;
    while (std::getline(initptclconf, initline))
    {
        initptclconf_list.push_back(std::stod(initline));
        // std::cout << initline << std::endl;
    }
    for(int num = 0; num < initptclconf_list.size()/4; num++)
    {
        ray_list.push_back(particle());
        ray_list.back().initptcl(initptclconf_list.at(0+(num*4)), initptclconf_list.at(1+(num*4)), initptclconf_list.at(2+(num*4)), initptclconf_list.at(3+(num*4)));
    }
    std::cout << "loaded " << ray_list.size() << " photon." << std::endl;
    while (std::getline(initscinticonf, initline))
    {
        initscinticonf_list.push_back(std::stod(initline));
        // std::cout << initline << std::endl;
    }
    for(int num = 0; num < initscinticonf_list.size()/9; num++)
    {
        scintillator.push_back(scinti());
        scintillator.back().initscinti(initscinticonf_list.at(0+(num*9)), initscinticonf_list.at(1+(num*9)), initscinticonf_list.at(2+(num*9)), initscinticonf_list.at(3+(num*9))*M_PI, initscinticonf_list.at(4+(num*9))*M_PI, initscinticonf_list.at(5+(num*9)), initscinticonf_list.at(6+(num*9)), initscinticonf_list.at(7+(num*9)), initscinticonf_list.at(8+(num*9)));
    }
    scintillator.back().initcs("./data/initcs_nai.conf");
    std::cout << "loaded " << scintillator.size() << " scintillator." << std::endl;
    std::cout << "please define run loop num" << std::endl;
    std::cin >> run_count;
    std::cout << "run loop defined " << run_count << std::endl;

    // 断面積の出力
    // crosssec_test(scintillator.back());

    // eigen並列処理無効化
    // Eigen::setNbThreads(1);
    // Eigen::initParallel();

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        double sum_ene = 0;
        std::cout << "----------\nrun: " << run << std::endl;
        bool break_flag = false;
        std::vector<particle> photon;
        photon.push_back(ray_list.at(0));
        photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
        // printf("thread = %d, run = %2d\n", omp_get_thread_num(), run);
        while (0 <= photon.back().ene_)
        {
            for (int scinti_num = 0; scinti_num < scintillator.size(); scinti_num++)
            {
                if (photon.back().ene_ <= 0)
                {
                    break_flag = true;
                    // break;
                }
                std::string outfilename = "./data/scinti_" + std::to_string(scinti_num) + ".dat";
                std::ofstream scinti_data(outfilename, std::ios::app);
                double traject_dist = scintillator.at(scinti_num).intersec_dist(ray_list.back(), photon.back()),
                    pe_cs = scintillator.at(scinti_num).crosssec(photon.back().ene_, 1),
                    pe_len = reactlen(pe_cs, scintillator.at(scinti_num).dens_),
                    cs_ang = cs_angle(photon.back().ene_),
                    cs_cs = scintillator.at(scinti_num).crosssec(photon.back().ene_, 2),
                    cs_len = reactlen(cs_cs, scintillator.at(scinti_num).dens_),
                    pp_cs = scintillator.at(scinti_num).crosssec(photon.back().ene_, 3),
                    pp_len = reactlen(pp_cs, scintillator.at(scinti_num).dens_);
                
                if(traject_dist > std::max({pe_len, cs_len, pp_len}) || traject_dist == -1)
                {
                    break_flag = true;
                    std::cout << "test1" << std::endl;
                }
                else
                {
                    // showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                    if(pe_len <= cs_len && pe_len <= pp_len)
                    {
                        sum_ene += photon.back().ene_;
                        break_flag = true;
                        std::cout << "test2" << std::endl;
                        // break;
                    }
                    else if(cs_len <= pe_len && cs_len <= pp_len)
                    {
                        sum_ene += photon.back().ene_ - scatphotonene(photon.back().ene_, cs_ang);
                        photon.back().ene_ = scatphotonene(photon.back().ene_, cs_ang);
                        photon.back().move(cs_len);
                        photon.back().turn(cs_ang);
                        std::cout << "test3" << std::endl;
                    }
                    else if (pp_len <= pe_len && pp_len <= cs_len)
                    {
                        photon.back().ene_ = 510.99895;
                        photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
                        photon.back().move(pp_len);
                        std::cout << "test4" << std::endl;
                    }
                    else
                    {
                        // std::cout << "error" << std::endl;
                        break_flag = true;
                        std::cout << "test5" << std::endl;
                        // break;
                    }
                    // std::cout << "sum_ene: " << sum_ene << std::endl;
                }

                // 2次反応以降
                if(photon.back().pt_x_ != ray_list.back().pt_x_ && photon.back().pt_y_ != ray_list.back().pt_y_ && photon.back().pt_z_ != ray_list.back().pt_z_)
                {
                    break_flag = true;
                }

                if (break_flag)
                {
                    if (sum_ene != 0)
                    {
                        std::cout << "sum_ene: " << sum_ene << std::endl;
                        #pragma omp critical
                        {
                            scinti_data << sum_ene << "\n";
                        }
                    }
                    std::cout << "test6" << std::endl;
                    sum_ene = 0;
                    // std::cout << "reset sum_ene: " << sum_ene << std::endl;
                    break;
                }
            }
            std::cout << "test7" << std::endl;
            if (break_flag)
            {
                break_flag = false;
                // std::cout << "test8" << std::endl;
                break;
            }
        }
    }
}