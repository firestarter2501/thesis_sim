#include "simfunc.h"

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

// void crosssec_test(scinti scintillator)
// {
//     std::ofstream cs_pe("../data/cs_pe.dat");
//     std::ofstream cs_cs("../data/cs_cs.dat");
//     std::ofstream cs_pp("../data/cs_pp.dat");

//     scintillator.initcs("../data/initcs_nai.conf");
// //std::cout << scintillator.crosssec_table_.size() << "\t" << scintillator.crosssec_table_.at(0).size() << std::endl;

//     // for(int i = 0; i < scintillator.crosssec_table_.size(); i++)
//     // {
//     //     for(int j = 0; j < scintillator.crosssec_table_.at(0).size(); j++)
//     //     {
//     //     std::cout << scintillator.crosssec_table_.at(i).at(j) << "\t";
//     //     }
//     // std::cout << "\n";
//     // }

// //std::cout << "------";

//     for(int ene = 10; ene < 90000; ene++)
//     {
//         cs_pe << ene << "\t" << scintillator.crosssec(ene, 1) << "\n";
//         cs_cs << ene << "\t" << scintillator.crosssec(ene, 2) << "\n";
//         cs_pp << ene << "\t" << scintillator.crosssec(ene, 3) << "\n";
//     }
// }

// Œõ“d‹zŽû’f–ÊÏ‚ð•Ô‚·ŠÖ”
double pe_crosssec(double ene, int z)
{
    double cs_tomscat = (8 * M_PI * pow(RELEC, 2)) / 3,
        fsconst = 1.0 / 137.0,
        crosssecperatom = 4 * pow(fsconst, 4) * sqrt(2) * pow(z, 5) * cs_tomscat * pow(pow(MEC2, 2) / ene, 7 / 2);
    return (pow(10, 14)/2.25)*crosssecperatom;
}

// ƒRƒ“ƒvƒgƒ“ŽU—‚·‚é‚Æ‚«‚ÌŽU—ŒõŽq‚ÌƒGƒlƒ‹ƒM[‚ð•Ô‚·ŠÖ”
double scatphotonene(double ene, double angle)
{
    double alpha = ene / MEC2;
    return ene / (1 + (alpha * (1 - std::cos(angle))));
}

// ƒNƒ‰ƒCƒ“m‰È‚ÌŽ®‚ðŒvŽZ‚·‚éŠÖ”
double kleinnishinaeq(double ene, double angle)
{
    double gamma = ene / MEC2;
    return (pow(10, 25)/3.5)*(pow(RELEC, 2) / 2) * (1 / pow(1 + gamma * (1 - std::cos(angle)), 2)) * (1 + pow(std::cos(angle), 2) + ((pow(gamma, 2) * pow(1 - std::cos(angle), 2)) / (1 + gamma * (1 - std::cos(angle)))));
    // return pow(RELEC, 2)*pow(1/(1+gamma*(1-std::cos(angle))), 2)*((1+pow(std::cos(angle), 2))/2)*(1+((pow(gamma, 2)*pow(1-std::cos(angle), 2))/((1+pow(std::cos(angle), 2))*(1+gamma*(1-std::cos(angle))))));
}

// von Neumann‚ÌŠü‹p–@‚ÅƒRƒ“ƒvƒgƒ“ŽU—‚ÌŠp“x•ª•z‚ÉŠî‚Ã‚¢‚½ŽU—Šp‚Ì—”‚ð•Ô‚·ŠÖ”
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
    // return initangle(randengine);
}

void cs_angle_test(int num)
{
    std::ofstream cs_angle_1kev_result("../data/cs_angle_1kev.dat"), cs_angle_100kev_result("../data/cs_angle_100kev.dat"), cs_angle_500kev_result("../data/cs_angle_500kev.dat"), cs_angle_2000kev_result("../data/cs_angle_2000kev.dat"), cs_angle_10000kev_result("../data/cs_angle_10000kev.dat");
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
    return -log(1 - initprob(randengine)) / (crosssec * dens * dens);
}

double normdist(double mean, double sdev)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::normal_distribution<> dist(mean, sdev);
    return std::abs(dist(engine));
}

void showinfo(std::vector<particle> photon, double traject_dist, double pe_len, double cs_len, double pp_len)
{
    std::cout << "-info-" << std::endl;
    std::cout << "photon_ene: " << photon.back().ene_ <<" photon_x: " << photon.back().pt_x_ << " photon_y: " << photon.back().pt_y_ << " photon_z: " << photon.back().pt_z_ <<  " photon_theta: " << photon.back().dir_theta_ << " photon_phi: " << photon.back().dir_phi_ << std::endl;
    std::cout << "traject_len: " << traject_dist << std::endl;
    std::cout << "pe_len: " << pe_len << " cs_len: " << cs_len<< " pp_len: " << pp_len << std::endl;
    std::cout << "---" << std::endl;
}

int randchoice(std::vector<double> probvec)
{
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::discrete_distribution<std::size_t> dist(
        probvec.begin(),
        probvec.end()
    );
    return dist(engine);
}

double lineareq(double ene, double slope, double intersec)
{
    return (ene*slope)+intersec;
}

void initcs(std::string conffilepath, std::vector<std::vector<double>> &crosssec_table)
{
    std::ifstream datafile(conffilepath);
    std::string line;
    while(getline(datafile, line))
    {
        std::istringstream stream(line);
        std::string field;
        std::vector<double> tmpvec;
        while(getline(stream, field, '\t'))
        {
            tmpvec.push_back(std::stod(field));
        }
        crosssec_table.push_back(tmpvec);
        tmpvec = {};
    }
}

double crosssec(double ene, int type, std::vector<std::vector<double>> crosssec_table)
{
    if (crosssec_table.at(crosssec_table.size()-1).at(0) <= ene)
    {
        return 0;
    }
    int ene_line = 0;
    while (crosssec_table.at(ene_line+1).at(0) < ene)
    {
        ene_line++;
    }
    
    // std::cout << "ene:" << ene << " ene_line:" << ene_line << std::endl;

    double alpha = (ene-crosssec_table.at(ene_line).at(0))/(crosssec_table.at(ene_line+1).at(0)-crosssec_table.at(ene_line).at(0)),
        cs_tmp = crosssec_table.at(ene_line).at(type)+((crosssec_table.at(ene_line+1).at(type)-crosssec_table.at(ene_line).at(type))*alpha);
    return cs_tmp;
}