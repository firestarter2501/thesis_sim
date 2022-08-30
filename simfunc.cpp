#include "simfunc.h"

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

// Œõ“d‹zû’f–ÊÏ‚ğ•Ô‚·ŠÖ”
double pe_crosssec(double ene, int z)
{
    double cs_tomscat = (8 * M_PI * pow(RELEC, 2)) / 3,
        fsconst = 1.0 / 137.0,
        crosssecperatom = 4 * pow(fsconst, 4) * sqrt(2) * pow(z, 5) * cs_tomscat * pow(pow(MEC2, 2) / ene, 7 / 2);
    return crosssecperatom;
}

// ƒRƒ“ƒvƒgƒ“U—‚·‚é‚Æ‚«‚ÌU—Œõq‚ÌƒGƒlƒ‹ƒM[‚ğ•Ô‚·ŠÖ”
double scatphotonene(double ene, double angle)
{
    double alpha = ene / MEC2;
    return ene / (1 + (alpha * (1 - std::cos(angle))));
}

// ƒNƒ‰ƒCƒ“m‰È‚Ì®‚ğŒvZ‚·‚éŠÖ”
double kleinnishinaeq(double ene, double angle)
{
    double gamma = ene / MEC2;
    return (pow(RELEC, 2) / 2) * (1 / pow(1 + gamma * (1 - std::cos(angle)), 2)) * (1 + pow(std::cos(angle), 2) + ((pow(gamma, 2) * pow(1 - std::cos(angle), 2)) / (1 + gamma * (1 - std::cos(angle)))));
    // return pow(RELEC, 2)*pow(1/(1+gamma*(1-std::cos(angle))), 2)*((1+pow(std::cos(angle), 2))/2)*(1+((pow(gamma, 2)*pow(1-std::cos(angle), 2))/((1+pow(std::cos(angle), 2))*(1+gamma*(1-std::cos(angle))))));
}

// von Neumann‚ÌŠü‹p–@‚ÅƒRƒ“ƒvƒgƒ“U—‚ÌŠp“x•ª•z‚ÉŠî‚Ã‚¢‚½U—Šp‚Ì—”‚ğ•Ô‚·ŠÖ”
double cs_angle(double ene)
{
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> initprob(0.0, 1.0);
    std::uniform_real_distribution<> initangle(0.0, M_PI);
    double max_crosssec = kleinnishinaeq(ene, 0.0),
        crosssec,
        randangle,
        randprob,
        prob;
    do
    {
        randangle = initangle(randengine);
        randprob = initprob(randengine);
        crosssec = kleinnishinaeq(ene, randangle);
        prob = 1.0 - (crosssec / max_crosssec);
    } while (randprob < prob);
    return randangle;
}

void cs_angle_test(int num)
{
    std::ofstream cs_angle_1kev_result("cs_angle_1kev.dat"), cs_angle_100kev_result("cs_angle_100kev.dat"), cs_angle_500kev_result("cs_angle_500kev.dat"), cs_angle_2000kev_result("cs_angle_2000kev.dat"), cs_angle_10000kev_result("cs_angle_10000kev.dat");
    std::vector<int> cs_angle_1kev, cs_angle_100kev, cs_angle_500kev, cs_angle_2000kev, cs_angle_10000kev;
    double dispparam = 10000;

    for (int arg = 0; arg <= 180; arg++)
    {
        cs_angle_1kev.push_back(0);
        cs_angle_100kev.push_back(0);
        cs_angle_500kev.push_back(0);
        cs_angle_2000kev.push_back(0);
        cs_angle_10000kev.push_back(0);
    }

    for (int i = 0; i < num; i++)
    {
        cs_angle_1kev.at(int(cs_angle(1) * (180.0 / M_PI))) += 1;
        cs_angle_100kev.at(int(cs_angle(100) * (180.0 / M_PI))) += 1;
        cs_angle_500kev.at(int(cs_angle(500) * (180.0 / M_PI))) += 1;
        cs_angle_2000kev.at(int(cs_angle(2000) * (180.0 / M_PI))) += 1;
        cs_angle_10000kev.at(int(cs_angle(10000) * (180.0 / M_PI))) += 1;
    }

    for (int arg = 0; arg < 180; arg++)
    {
        cs_angle_1kev_result << arg + 1 << "\t" << dispparam * cs_angle_1kev.at(arg) / cs_angle_1kev.at(0) << "\n";
        cs_angle_100kev_result << arg + 1 << "\t" << dispparam * cs_angle_100kev.at(arg) / cs_angle_100kev.at(0) << "\n";
        cs_angle_500kev_result << arg + 1 << "\t" << dispparam * cs_angle_500kev.at(arg) / cs_angle_500kev.at(0) << "\n";
        cs_angle_2000kev_result << arg + 1 << "\t" << dispparam * cs_angle_2000kev.at(arg) / cs_angle_2000kev.at(0) << "\n";
        cs_angle_10000kev_result << arg + 1 << "\t" << dispparam * cs_angle_10000kev.at(arg) / cs_angle_10000kev.at(0) << "\n";
    }
}

double reactlen(double crosssec, double dens)
{
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> initprob(0.0, 1.0);
    return -log(1 - initprob(randengine)) / (crosssec * dens);
}