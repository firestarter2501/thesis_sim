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
    // strbuf = std::cout.rdbuf(ofstr.rdbuf());
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
    std::cin >> run_count;

    /*出力ファイル定義*/
    std::ofstream scintiofstr;
    scintiofstr.open("../data/scinti.dat", std::ios::app);

    #pragma omp parallel for
    for (int run = 0; run < run_count; run++)
    {
        particle ptcl = ray_list.back();
        ptcl.initptcl(ptcl.ene_, ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_);

        /*オブジェクトとロジックを配置*/
        scinti scinti1;
        scinti1.initscinti("scinti1", "initcs_nai");
        double scinti1tmp = 0;

        while (ptcl.ene_ > 0)
        {
            scinti1tmp = scinti1.scintillation(ptcl);
            
            if (scinti1tmp > 0)
            {
                #pragma omp critical
                {
                scintiofstr << scinti1tmp << "\n";
                }
            }
        }
    }
    scintiofstr.close();
}