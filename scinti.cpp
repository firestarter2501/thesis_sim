#include "scinti.h"

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23) // mol^-1

std::vector<bool> scinti::intersecenablecheck(std::vector<std::vector<double>> intersecvec, particle ptcl)
{
    std::vector<bool> returnvec;
    for (int i = 0; i < intersecvec.size(); i++)
    {
        Eigen::Vector3d intersecpoint, scintipoint, centerdist, ptcldir, intersecdir, dirnorm;
        intersecpoint << intersecvec.at(i).at(0), intersecvec.at(i).at(1), intersecvec.at(i).at(2);
        scintipoint << this->pt_x_, this->pt_y_, this->pt_z_;
        centerdist = intersecpoint - scintipoint;

        ptcldir << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
        intersecdir << intersecvec.at(i).at(0) - ptcl.pt_x_, intersecvec.at(i).at(1) - ptcl.pt_y_, intersecvec.at(i).at(2) - ptcl.pt_z_;
        ptcldir.normalize();
        intersecdir.normalize();
        dirnorm = intersecdir - ptcldir;
        std::cout << "dirnorm.norm(): " << dirnorm.norm() << std::endl;
        if (std::abs(centerdist.norm()) <= std::sqrt(2)*this->rad_ && dirnorm.norm() < 0.000001)
        {
            returnvec.push_back(true);
        }
        else
        {
            returnvec.push_back(false);
        }
    }
    return returnvec;
}

double scinti::ptclinsidecheck(particle ptcl)
{
    Eigen::Vector3d ptclpoint, scintipoint, centerdist;
    ptclpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
    scintipoint << this->pt_x_, this->pt_y_, this->pt_z_;
    centerdist = ptclpoint - scintipoint;

    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<bool> enablecheck = intersecenablecheck(intersec, ptcl);
    int ec_truecount = std::count(enablecheck.begin(), enablecheck.end(), true);

    if (std::abs(centerdist.norm()) <= std::sqrt(2)*this->rad_ && ec_truecount != 2)
    {
        return 0;
    }
    else
    {
        auto ec_true1 = std::find(enablecheck.begin(), enablecheck.end(), true);
        int type1 = std::distance(enablecheck.begin(), ec_true1);
        auto ec_true2 = std::find(++ec_true1, enablecheck.end(), true);
        int type2 = std::distance(enablecheck.begin(), ec_true2);
        Eigen::Vector3d type1vec, type2vec, traject1vec, traject2vec;
        type1vec << intersec.at(type1).at(0), intersec.at(type1).at(1), intersec.at(type1).at(2);
        type2vec << intersec.at(type2).at(0), intersec.at(type2).at(1), intersec.at(type2).at(2);
        traject1vec = type1vec - ptclpoint;
        traject2vec = type2vec - ptclpoint;
        std::cout << "ptclinsidecheck: " << std::min({ std::abs(traject1vec.norm()), std::abs(traject2vec.norm()) }) << std::endl;
        return std::min({ std::abs(traject1vec.norm()), std::abs(traject2vec.norm()) });
    }
}

std::string scinti::showfacetype(int type)
{
    if (type == 0)
    {
        return "front";
    }
    else if (type == 1)
    {
        return "back";
    }
    else if (type == 2)
    {
        return "side1";
    }
    else if (type == 3)
    {
        return "side2";
    }
    else
    {
        return "facetypeerror";
    }
}

double scinti::intersec_dist(particle ptcl)
{
    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<bool> enablecheck = intersecenablecheck(intersec, ptcl);

    // intersecの座標表示
    for (int i = 0; i < intersec.size(); i++)
    {
        std::cout << showfacetype(i) << ": ";
        for (int j = 0; j < intersec.at(0).size(); j++)
        {
            std::cout << intersec.at(i).at(j) << "\t";
        }
        std::cout << enablecheck.at(i) << "\n";
    }

    int ec_truecount = std::count(enablecheck.begin(), enablecheck.end(), true);
    if (ec_truecount == 0)
    {
        std::cout << "outside(all false)" << std::endl;
        return -1;
    }
    else if (ec_truecount == 1)
    {
        auto ec_true = std::find(enablecheck.begin(), enablecheck.end(), true);
        int type = std::distance(enablecheck.begin(), ec_true);
        Eigen::Vector3d ptclvec, intersecvec, trajectvec;
        ptclvec << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
        intersecvec << intersec.at(type).at(0), intersec.at(type).at(1), intersec.at(type).at(2);
        trajectvec = intersecvec - ptclvec;
        std::cout << "ptcl & " << showfacetype(type) << " dist: " << trajectvec.norm() << std::endl;
        return trajectvec.norm();
    }
    else if (ec_truecount == 2)
    {
        auto ec_true1 = std::find(enablecheck.begin(), enablecheck.end(), true);
        int type1 = std::distance(enablecheck.begin(), ec_true1);
        auto ec_true2 = std::find(++ec_true1, enablecheck.end(), true);
        int type2 = std::distance(enablecheck.begin(), ec_true2);
        Eigen::Vector3d type1vec, type2vec, trajectvec;
        type1vec << intersec.at(type1).at(0), intersec.at(type1).at(1), intersec.at(type1).at(2);
        type2vec << intersec.at(type2).at(0), intersec.at(type2).at(1), intersec.at(type2).at(2);
        trajectvec = type2vec - type1vec;
        std::cout << showfacetype(type1) << " & " << showfacetype(type2) << " dist: " << trajectvec.norm() << std::endl;
        return trajectvec.norm();
    }
    else
    {
        std::cout << "error(too many true)" << std::endl;
        return -1;
    }
}

void scinti::initcs(std::string conffilepath)
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

double scinti::crosssec(double ene, int type)
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

void scinti::initscinti(double pt_x, double pt_y, double pt_z, double theta, double phi, double depth, double z, double dens, double atomweight)
{
    this->pt_x_ = pt_x;
    this->pt_y_ = pt_y;
    this->pt_z_ = pt_z;
    this->dir_theta_ = theta;
    this->dir_phi_ = phi;
    this->depth_ = depth;
    this->rad_ = depth / 2;
    this->z_ = z;
    this->dens_ = dens;
    this->ndens_ = (dens / atomweight) * NUMA;
    this->atomweight_ = atomweight;
}

std::vector<std::vector<double>> scinti::intersec(particle ptcl)
{
    std::vector<std::vector<double>> return_point = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
    Eigen::Vector3d traject, ptclpoint, scinti_centerline, scinti_frontcenter, scinti_front_intersec, frontdist, scinti_backcenter, scinti_back_intersec, backdist;

    traject << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
    traject.normalize();

    ptclpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;

    // std::cout << "dir_theta_: " << this->dir_theta_ << ", dir_phi_: " << this->dir_phi_ << std::endl;

    scinti_centerline << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_)* std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    scinti_centerline.normalize();
    // std::cout << "scinti_ceneterline: " << scinti_centerline(0) << ", " << scinti_centerline(1) << ", " << scinti_centerline(2) << std::endl;

    scinti_frontcenter << this->pt_x_ + (this->depth_/2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (this->depth_ / 2) * std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (this->depth_ / 2) * std::cos(this->dir_theta_);
    // std::cout << "scinti_frontcenter: " << scinti_frontcenter(0) << ", " << scinti_frontcenter(1) << ", " << scinti_frontcenter(2) << std::endl;

    double front_t = (scinti_centerline(0)*(scinti_frontcenter(0)-ptcl.pt_x_) + scinti_centerline(1) * (scinti_frontcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_frontcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    scinti_front_intersec << ptcl.pt_x_ + traject(0) * front_t, ptcl.pt_y_ + traject(1) * front_t, ptcl.pt_z_ + traject(2) * front_t;

    frontdist = scinti_front_intersec - scinti_frontcenter;
    if (std::abs(frontdist.norm()) <= this->rad_)
    {
        return_point.at(0).at(0) = scinti_front_intersec(0);
        return_point.at(0).at(1) = scinti_front_intersec(1);
        return_point.at(0).at(2) = scinti_front_intersec(2);
        // std::cout << "front: " << return_point.at(0).at(0) << ", " << return_point.at(0).at(1) << ", " << return_point.at(0).at(2) << std::endl;
    }

    scinti_backcenter << this->pt_x_ + (-this->depth_ / 2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (-this->depth_ / 2)* std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (-this->depth_ / 2)* std::cos(this->dir_theta_);
    // std::cout << "scinti_backcenter: " << scinti_backcenter(0) << ", " << scinti_backcenter(1) << ", " << scinti_backcenter(2) << std::endl;

    double back_t = (scinti_centerline(0) * (scinti_backcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_backcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_backcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    scinti_back_intersec << ptcl.pt_x_ + traject(0) * back_t, ptcl.pt_y_ + traject(1) * back_t, ptcl.pt_z_ + traject(2) * back_t;

    backdist = scinti_back_intersec - scinti_backcenter;
    if (std::abs(backdist.norm()) <= this->rad_)
    {
        return_point.at(1).at(0) = scinti_back_intersec(0);
        return_point.at(1).at(1) = scinti_back_intersec(1);
        return_point.at(1).at(2) = scinti_back_intersec(2);
        // std::cout << "back: " << return_point.at(1).at(0) << ", " << return_point.at(1).at(1) << ", " << return_point.at(1).at(2) << std::endl;
    }

    Eigen::Vector3d p, p2, vecs, v;   
    p << scinti_frontcenter(0) - ptcl.pt_x_, scinti_frontcenter(1) - ptcl.pt_y_, scinti_frontcenter(2) - ptcl.pt_z_;
    p2 << scinti_backcenter(0) - ptcl.pt_x_, scinti_backcenter(1) - ptcl.pt_y_, scinti_backcenter(2) - ptcl.pt_z_;
    vecs << p2(0) - p(0), p2(1) - p(1), p2(2) - p(2);
    v = traject;
    double Dvv = v.dot(v),
        Dsv = vecs.dot(v),
        Dpv = p.dot(v),
        Dss = vecs.dot(vecs),
        Dps = p.dot(vecs),
        Dpp = p.dot(p),
        rr = this->rad_*this->rad_;
    if (Dss == 0)
    {
        // std::cout << "void cilinder" << std::endl;
        return return_point;
    }
    double A = Dvv - Dsv * Dsv / Dss,
        B = Dpv - Dps * Dsv / Dss,
        C = Dpp - Dps * Dps / Dss - rr; 
    if (A == 0)
    {
        return return_point;
    }
    double ds = B * B - A * C;
    if (ds < 0)
    {
        // std::cout << "side no collision" << std::endl;
        return return_point;
    }
    ds = sqrt(ds);
    double a1 = (B - ds) / A,
        a2 = (B + ds) / A;
    return_point.at(2).at(0) = ptcl.pt_x_ + a1 * v(0);
    return_point.at(2).at(1) = ptcl.pt_y_ + a1 * v(1);
    return_point.at(2).at(2) = ptcl.pt_z_ + a1 * v(2);
    return_point.at(3).at(0) = ptcl.pt_x_ + a2 * v(0);
    return_point.at(3).at(1) = ptcl.pt_y_ + a2 * v(1);
    return_point.at(3).at(2) = ptcl.pt_z_ + a2 * v(2);

    // std::cout << "side1: " << return_point.at(2).at(0) << ", " << return_point.at(2).at(1) << ", " << return_point.at(2).at(2) << std::endl;
    // std::cout << "side2: " << return_point.at(3).at(0) << ", " << return_point.at(3).at(1) << ", " << return_point.at(3).at(2) << std::endl;
    
    return return_point;
}