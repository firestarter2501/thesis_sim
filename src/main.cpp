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
    std::ofstream scinti1ofstr("../data/scinti1.dat");
    std::ofstream scinti2ofstr("../data/scinti2.dat");

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        //*std::cout << std::endl << "-----------------------\nrun: " << run << std::endl;
        particle ptcl = ray_list.back();
        ptcl.initptcl(ptcl.ene_, ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_);

        /*オブジェクト初期化*/
        scinti scinti1;
        double scinti1tmp = 1;

        scinti scinti2;
        double scinti2tmp = 1;

        while (scinti1tmp+scinti2tmp > -2)
        {
            if (ptcl.ene_ <= 0)
            {
                break;
            }

            /*オブジェクトとロジックを配置*/
            //*std::cout << "-----scinti1-----" << std::endl;
            scinti1tmp = scinti1.scintillation("scinti1", "initcs_nai", ptcl);

            if (ptcl.ene_ <= 0)
            {
                break;
            }

            //*std::cout << "-----scinti2-----" << std::endl;
            scinti2tmp = scinti2.scintillation("scinti2", "initcs_nai", ptcl);

            //*std::cout << "scinti1tmp: " << scinti1tmp << ", scinti2tmp: " << scinti2tmp << std::endl;

            if (scinti1tmp > 0 && scinti2tmp > 0)
            {
                //*std::cout << "enable" << std::endl;
                #pragma omp critical
                {
                scinti1ofstr << scinti1tmp << "\n";
                scinti2ofstr << scinti2tmp << "\n";
                }
            }
        }
    }
}