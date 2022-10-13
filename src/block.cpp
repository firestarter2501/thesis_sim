#include "block.h"

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23) // mol^-1

using vsize_t = std::vector<int>::size_type;

void block::initcs(std::string conffilepath)
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
        this->crosssec_table_.push_back(tmpvec);
        tmpvec = {};
    }
}

double block::crosssec(double ene, int type)
{
    if (this->crosssec_table_.at(this->crosssec_table_.size()-1).at(0) <= ene)
    {
        return 0;
    }
    int ene_line = 0;
    while (this->crosssec_table_.at(ene_line+1).at(0) < ene)
    {
        ene_line++;
    }
    
    // std::cout << "ene:" << ene << " ene_line:" << ene_line << std::endl;

    double alpha = (ene-this->crosssec_table_.at(ene_line).at(0))/(this->crosssec_table_.at(ene_line+1).at(0)-this->crosssec_table_.at(ene_line).at(0)),
        cs_tmp = this->crosssec_table_.at(ene_line).at(type)+((this->crosssec_table_.at(ene_line+1).at(type)-this->crosssec_table_.at(ene_line).at(type))*alpha);
    return cs_tmp;

}