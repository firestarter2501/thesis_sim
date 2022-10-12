#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include "simfunc.h"
#include "particle.h"
#include "scinti.h"

using vsize_t = std::vector<int>::size_type;

int main()
{
    std::vector<particle> ray_list;
    std::vector<double> initptclconf_list;
    std::ifstream initptclconf("./data/initptcl.conf");
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
    std::cout << "please define run loop num" << std::endl;
    std::cin >> run_count;
    std::cout << "run loop defined " << run_count << std::endl;

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        std::cout << std::endl << "-----------------------\nrun: " << run << std::endl;
        particle ptcl = ray_list.back();
        ptcl.initptcl(ptcl.ene_, ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_);
        std::vector<bool> react_flag;
        bool tmp_react_flag = true;
        bool absorp_flag = false;
        /*以下に使用するオブジェクトを定義*/
        scinti scinti1;
        react_flag.push_back(true);

        scinti scinti2;
        react_flag.push_back(true);

        /*---------------------------*/

        while (0 < std::count(react_flag.begin(), react_flag.end(), true))
        {
            if (absorp_flag)
            {
                break;
            }
            /*以下に使用するオブジェクトを配置*/
            scinti1.scintillation("scinti1", "initcs_nai", ptcl, tmp_react_flag, absorp_flag);
            react_flag.at(0) = tmp_react_flag;

            scinti1.scintillation("scinti2", "initcs_nai", ptcl, tmp_react_flag, absorp_flag);
            react_flag.at(1) = tmp_react_flag;

            /*---------------------------*/
        }
    }
}