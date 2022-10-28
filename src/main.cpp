#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include "simfunc.h"
#include "particle.h"
#include "scinti.h"
#include "block.h"

using vsize_t = std::vector<int>::size_type;

int main()
{
    std::vector<particle> ray_list;
    std::vector<double> initptclconf_list;
    std::ifstream initptclconf("../data/initptcl.conf");
    std::ofstream ofstr("../data/sim.log");
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
    //*std::cout << "loaded " << ray_list.size() << " photon." << std::endl;
    //*std::cout << "please define run loop num" << std::endl;
    std::cin >> run_count;
    //*std::cout << "run loop defined " << run_count << std::endl;

    /*出力ファイル定義*/
    std::ofstream scintiofstr;
    scintiofstr.open("../data/scinti.dat", std::ios::app);
    // std::ofstream scinti1ofstr("../data/scinti1.dat");
    // std::ofstream scinti2ofstr("../data/scinti2.dat");

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        //*std::cout << std::endl << "-----------------------\nrun: " << run << std::endl;
        particle ptcl = ray_list.back();
        ptcl.initptcl(ptcl.ene_, ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_);

        /*オブジェクトとロジックを配置*/
        scinti scinti1;
        scinti1.initscinti("scinti1", "initcs_nai");
        double scinti1_dist = 0, scinti1tmp = 0;

        scinti scinti2;
        scinti2.initscinti("scinti2_150deg", "initcs_nai");
        double scinti2_dist = 0, scinti2tmp = 0;

        while (scinti1_dist + scinti2_dist > -2 && ptcl.ene_ > 0)
        {
            scinti1_dist = scinti1.intersec_dist(ptcl);
            scinti2_dist = scinti2.intersec_dist(ptcl);
            //*std::cout << "scinti1.intersec_dist(ptcl): " << scinti1_dist << ", scinti2.intersec_dist(ptcl): " << scinti2_dist << std::endl;
            if (0 < scinti1_dist && scinti2_dist <= 0)
            {
                //*std::cout << "-----scinti1-----" << std::endl;
                scinti1tmp = scinti1.scintillation(ptcl);

                //*std::cout << "-----scinti2-----" << std::endl;
                scinti2tmp = scinti2.scintillation(ptcl);
            }
            else if (scinti1_dist <= 0 && 0 < scinti2_dist)
            {
                //*std::cout << "-----scinti2-----" << std::endl;
                scinti2tmp = scinti2.scintillation(ptcl);

                //*std::cout << "-----scinti1-----" << std::endl;
                scinti1tmp = scinti1.scintillation(ptcl);
            }
            else
            {
                //*std::cout << "error" << std::endl;
                break;
            }
            
            if (scinti1tmp > 0 && scinti2tmp > 0)
            {
                //*std::cout << "enable data" << std::endl;
                #pragma omp critical
                {
                scintiofstr << scinti1tmp << "\t" << scinti2tmp << "\n";
                // scinti1ofstr << scinti1tmp << "\n";
                // scinti2ofstr << scinti2tmp << "\n";
                }
            }
        }
    }
    scintiofstr.close();
}