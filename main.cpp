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
#define PMTSDEVSLOPE 0.018471
#define PMTSDEVINTERSEC 12.2998905

using vsize_t = std::vector<int>::size_type;

int main()
{
    std::vector<particle> ray_list;
    std::vector<scinti> scintillator;
    std::vector<double> initptclconf_list;
    std::vector<double> initscinticonf_list;
    std::ifstream initptclconf("./data/initptcl.conf");
    std::ifstream initscinticonf("./data/initscinti.conf");
    // std::ifstream initscinticonf("./data/initscinti2.conf");
    std::ofstream ofstr("./data/sim.log");
    // std::streambuf* strbuf;
    std::string initline;
    /*strbuf = */std::cout.rdbuf(ofstr.rdbuf());
    int run_count = 0;
    while (std::getline(initptclconf, initline))
    {
        initptclconf_list.push_back(std::stod(initline));
        // std::cout << initline << std::endl;
    }
    for(vsize_t num = 0; num < initptclconf_list.size()/4; num++)
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
    for(vsize_t num = 0; num < initscinticonf_list.size()/9; num++)
    {
        scintillator.push_back(scinti());
        scintillator.back().initcs("./data/initcs_nai.conf");
        scintillator.back().initscinti(initscinticonf_list.at(0+(num*9)), initscinticonf_list.at(1+(num*9)), initscinticonf_list.at(2+(num*9)), initscinticonf_list.at(3+(num*9))*M_PI, initscinticonf_list.at(4+(num*9))*M_PI, initscinticonf_list.at(5+(num*9)), initscinticonf_list.at(6+(num*9)), initscinticonf_list.at(7+(num*9)), initscinticonf_list.at(8+(num*9)));
    }
    std::cout << "loaded " << scintillator.size() << " scintillator." << std::endl;
    std::cout << "please define run loop num" << std::endl;
    std::cin >> run_count;
    std::cout << "run loop defined " << run_count << std::endl;

    // 断面積の出力
    crosssec_test(scintillator.back());
    cs_angle_test(10000);

    // eigen並列処理無効化
    // Eigen::setNbThreads(1);
    // Eigen::initParallel();

    particle test;
    test.initptcl(661.6, 0, 0, 0);
    test.turn_test(M_PI*1/3, M_PI*3/4, M_PI*3/6);

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        int reactcount = 0;
        std::cout << std::endl << "-----------------------\nrun: " << run << std::endl;
        std::vector<bool> scinti_react(scintillator.size(), true);
        std::vector<scinti> scinti_inloop = scintillator;
        std::vector<particle> photon;
        photon.push_back(ray_list.at(0));
        photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
        // printf("thread = %d, run = %2d\n", omp_get_thread_num(), run);
        while (0 < std::count(scinti_react.begin(), scinti_react.end(), true))
        {
            for (vsize_t scinti_num = 0; scinti_num < scinti_inloop.size(); scinti_num++)
            {
                std::cout << std::endl << "--------------\nscinti_num: " << scinti_num << std::endl;
                std::string outfilename = "./data/scinti_" + std::to_string(scinti_num) + ".dat";
                std::ofstream scinti_data(outfilename, std::ios::app);
                std::cout << "beforetraject_reactcount: " << reactcount << std::endl;
                double traject_dist = scinti_inloop.at(scinti_num).intersec_dist(photon.back()),
                    pe_cs = scinti_inloop.at(scinti_num).crosssec(photon.back().ene_, 1),
                    pe_len = reactlen(pe_cs, scinti_inloop.at(scinti_num).dens_),
                    cs_ang = cs_angle(photon.back().ene_),
                    cs_cs = scinti_inloop.at(scinti_num).crosssec(photon.back().ene_, 2),
                    cs_len = reactlen(cs_cs, scinti_inloop.at(scinti_num).dens_),
                    pp_cs = scinti_inloop.at(scinti_num).crosssec(photon.back().ene_, 3),
                    pp_len = reactlen(pp_cs, scinti_inloop.at(scinti_num).dens_);
                
                if(traject_dist < std::min({pe_len, cs_len, pp_len}) || traject_dist == -1/* || (scinti_num == 1 && photon.back().ene_ == ray_list.back().ene_)*/)
                {
                    scinti_react.at(scinti_num) = false;
                    std::cout << "outside or too short(traject_dist: " << traject_dist << ", pe_len: " << pe_len << ", cs_len: " << cs_len << ", pp_len: " << pp_len << std::endl;
                    showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                    std::cout << "ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << std::endl;
                }
                else
                {
                    std::cout << "before ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << std::endl;
                    if(pe_len <= cs_len && pe_len <= pp_len)
                    {
                        scinti_inloop.at(scinti_num).ene_buffer_ += normdist(photon.back().ene_, lineareq(photon.back().ene_, PMTSDEVSLOPE, PMTSDEVINTERSEC));
                        reactcount++;
                        scinti_react.at(scinti_num) = false;
                        std::cout << "---pe---" << std::endl;
                        std::cout << "ene_buffer: " << scinti_inloop.at(scinti_num).ene_buffer_ << std::endl;
                        showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                    }
                    else if(cs_len <= pe_len && cs_len <= pp_len)
                    {
                        showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                        Eigen::Vector3d beforepoint, afterpoint, moveddist;
                        scinti_inloop.at(scinti_num).ene_buffer_ += normdist(photon.back().ene_ - scatphotonene(photon.back().ene_, cs_ang), lineareq(photon.back().ene_ - scatphotonene(photon.back().ene_, cs_ang), PMTSDEVSLOPE, PMTSDEVINTERSEC));
                        photon.back().ene_ = scatphotonene(photon.back().ene_, cs_ang);
                        beforepoint << photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_;
                        photon.back().move(scinti_inloop.at(scinti_num).ptclinsidecheck(photon.back()) + cs_len);
                        afterpoint << photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_;
                        moveddist = afterpoint - beforepoint;
                        std::cout << "moved dist: " << std::abs(moveddist.norm()) << std::endl;
                        photon.back().turn(cs_ang);
                        reactcount++;
                        scinti_react.at(scinti_num) = true;
                        std::cout << "---cs---" << std::endl;
                        std::cout << "ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << " cs_len: " << cs_len << " cs_ang: " << cs_ang  << std::endl;
                        showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                    }
                    else if (pp_len <= pe_len && pp_len <= cs_len)
                    {
                        Eigen::Vector3d beforepoint, afterpoint, moveddist;
                        photon.back().ene_ = MEC2;
                        photon.back().initptcl(photon.back().ene_, photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_);
                        beforepoint << photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_;
                        photon.back().move(scinti_inloop.at(scinti_num).ptclinsidecheck(photon.back()) + pp_len);
                        afterpoint << photon.back().pt_x_, photon.back().pt_y_, photon.back().pt_z_;
                        moveddist = afterpoint - beforepoint;
                        std::cout << "moved dist: " << std::abs(moveddist.norm()) << std::endl;
                        // reactcount++;
                        scinti_react.at(scinti_num) = true;
                        std::cout << "---pp---" << std::endl;
                        std::cout << "ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << " pp_len: " << pp_len << std::endl;
                        showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                    }
                    else
                    {
                        std::cout << "judge error" << std::endl;
                        showinfo(photon, traject_dist, pe_len, cs_len, pp_len);
                        scinti_react.at(scinti_num) = false;
                    }
                }

                // 2次反応以降無効化
                // if(reactcount >= 1)
                // {
                //     break_flag = true;
                //     // std::cout << "break for first" << std::endl;
                // }

                if (!scinti_react.at(scinti_num))
                {
                    #pragma omp critical
                    {
                    if (0 < scinti_inloop.at(scinti_num).ene_buffer_)
                    {
                        std::cout << "final ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << " reactcount: " << reactcount << std::endl;
                        scinti_data << scinti_inloop.at(scinti_num).ene_buffer_ << "\n";
                        std::cout << "sum_ene check&add done" << std::endl;
                    }
                    }
                    scinti_inloop.at(scinti_num).ene_buffer_ = 0;
                    reactcount = 0;
                    std::cout << "reset ene_buffer_: " << scinti_inloop.at(scinti_num).ene_buffer_ << " reactcount: " << reactcount << std::endl;
                    // break;
                }
            }
            std::cout << "scinti for break" << std::endl;
            // if (std::count(scinti_react.begin(), scinti_react.end(), true) == 0)
            // {
            //     std::cout << "while break" << std::endl;
            //     break;
            // }
        }
    }
}