#include "scinti.h"

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23) // mol^-1

double scinti::limtozero(double num)
{
    if (0 < num && num < 0.00000001)
    {
        return 0;
    }
    else
    {
        return num;
    }
}

std::vector<double> scinti::scintibc()
{
    double f_l_x = this->pt_x_ - (sqrt(2) * this->rad_ * sin((M_PI/4)-(M_PI-this->dir_phi_))),
        f_r_x = this->pt_x_ + (sqrt(2) * this->rad_ * sin((M_PI/4)+(M_PI-this->dir_phi_))),
        f_l_y = this->pt_y_ - (sqrt(2) * this->rad_ * sin((M_PI/2)-(M_PI-this->dir_phi_)-(M_PI/4))),
        f_r_y = this->pt_y_ - (sqrt(2) * this->rad_ * sin((M_PI-this->dir_phi_)-(M_PI/4))),
        f_t_z = this->pt_z_ + (sqrt(2) * this->rad_ * sin((M_PI/4)-((M_PI/2)-this->dir_theta_))),
        f_u_z = this->pt_z_ - (sqrt(2) * this->rad_ * sin((M_PI/4)+((M_PI/2)-this->dir_theta_))),
        b_l_x = this->pt_x_ - (sqrt(2) * this->rad_ * sin((M_PI/4)+(M_PI-this->dir_phi_))),
        b_r_x = this->pt_x_ + (sqrt(2) * this->rad_ * sin((M_PI/4)-(M_PI-this->dir_phi_))),
        b_l_y = this->pt_y_ + (sqrt(2) * this->rad_ * sin((M_PI-this->dir_phi_)-(M_PI/4))),
        b_r_y = this->pt_y_ + (sqrt(2) * this->rad_ * sin((M_PI/2)-(M_PI-this->dir_phi_)-(M_PI/4))),
        b_t_z = this->pt_z_ + (sqrt(2) * this->rad_ * sin((M_PI/4)+((M_PI/2)-this->dir_theta_))),
        b_u_z = this->pt_z_ - (sqrt(2) * this->rad_ * sin((M_PI/4)-((M_PI/2)-this->dir_theta_)));
    std::vector<double> return_vec = {f_l_x, f_r_x, f_l_y, f_r_y, f_t_z, f_u_z, b_l_x, b_r_x, b_l_y, b_r_y, b_t_z, b_u_z};
    return return_vec;
}

bool scinti::initpointcheck(particle initptcl, particle ptcl)
{
    if (ptcl.pt_x_ == initptcl.pt_x_ && ptcl.pt_y_ == initptcl.pt_y_ && ptcl.pt_z_ == initptcl.pt_z_)
    {
        return true;
    }
    else{
        return false;
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
    // std::cout << "test1" << std::endl;
    if(this->crosssec_table_.at(this->crosssec_table_.size()-1).at(0) <= ene)
    {
        return 0;
    }
    int ene_line = 0;
    // std::cout << "test2" << std::endl;
    while(this->crosssec_table_.at(ene_line).at(0) < ene)
    {
        // std::cout << "test3" << std::endl;
        ene_line++;
    }

    // std::cout << "test4" << std::endl;

    double cs_tmp = this->crosssec_table_.at(ene_line).at(type)+((this->crosssec_table_.at(ene_line+1).at(type)-this->crosssec_table_.at(ene_line).at(type))*(ene-this->crosssec_table_.at(ene_line).at(0))/(this->crosssec_table_.at(ene_line+1).at(0)-this->crosssec_table_.at(ene_line).at(0)));

    if(cs_tmp < 0.001)
    {
        return 0;
    }
    else
    {
        return cs_tmp;
    }
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
    Eigen::Vector3d traject, scinti_centerline, scinti_frontcenter, scinti_front_intersec, scinti_backcenter, scinti_back_intersec;

    traject << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
    traject.normalize();

    // std::cout << "dir_theta_: " << this->dir_theta_ << ", dir_phi_: " << this->dir_phi_ << std::endl;

    scinti_centerline << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_)* std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    scinti_centerline.normalize();
    // std::cout << "scinti_ceneterline: " << scinti_centerline(0) << ", " << scinti_centerline(1) << ", " << scinti_centerline(2) << std::endl;

    scinti_frontcenter << this->pt_x_ + (this->depth_/2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (this->depth_ / 2) * std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (this->depth_ / 2) * std::cos(this->dir_theta_);
    // std::cout << "scinti_frontcenter: " << scinti_frontcenter(0) << ", " << scinti_frontcenter(1) << ", " << scinti_frontcenter(2) << std::endl;

    double front_t = (scinti_centerline(0)*(scinti_frontcenter(0)-ptcl.pt_x_) + scinti_centerline(1) * (scinti_frontcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_frontcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (std::sqrt(std::pow(scinti_frontcenter(0) - (ptcl.pt_x_ + traject(0)*front_t), 2) + std::pow(scinti_frontcenter(1) - (ptcl.pt_y_ + traject(1) * front_t), 2) + std::pow(scinti_frontcenter(2) - (ptcl.pt_z_ + traject(2) * front_t), 2)) < this->rad_)
    {
        double fx = limtozero(ptcl.pt_x_ + traject(0) * front_t),
            fy = limtozero(ptcl.pt_y_ + traject(1) * front_t),
            fz = limtozero(ptcl.pt_z_ + traject(2) * front_t);
        return_point.at(0).at(0) = fx;
        return_point.at(0).at(1) = fy;
        return_point.at(0).at(2) = fz;
        std::cout << "front: " << return_point.at(0).at(0) << ", " << return_point.at(0).at(1) << ", " << return_point.at(0).at(2) << std::endl;
    }

    scinti_backcenter << this->pt_x_ + (-this->depth_ / 2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (-this->depth_ / 2)* std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (-this->depth_ / 2)* std::cos(this->dir_theta_);
    // std::cout << "scinti_backcenter: " << scinti_backcenter(0) << ", " << scinti_backcenter(1) << ", " << scinti_backcenter(2) << std::endl;
    double back_t = (scinti_centerline(0) * (scinti_backcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_backcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_backcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (std::sqrt(std::pow(scinti_backcenter(0) - (ptcl.pt_x_ + traject(0) * back_t), 2) + std::pow(scinti_backcenter(1) - (ptcl.pt_y_ + traject(1) * back_t), 2) + std::pow(scinti_backcenter(2) - (ptcl.pt_z_ + traject(2) * back_t), 2)) < this->rad_)
    {
        double bx = limtozero(ptcl.pt_x_ + traject(0) * back_t),
            by = limtozero(ptcl.pt_y_ + traject(1) * back_t),
            bz = limtozero(ptcl.pt_z_ + traject(2) * back_t);
        return_point.at(1).at(0) = limtozero(ptcl.pt_x_ + traject(0) * back_t);
        return_point.at(1).at(1) = limtozero(ptcl.pt_y_ + traject(1) * back_t);
        return_point.at(1).at(2) = limtozero(ptcl.pt_z_ + traject(2) * back_t);
        std::cout << "back: " << return_point.at(1).at(0) << ", " << return_point.at(1).at(1) << ", " << return_point.at(1).at(2) << std::endl;
    }

    Eigen::Vector3d p, p2, vecs, v;   
    p << scinti_front_intersec(0) - ptcl.pt_x_, scinti_front_intersec(1) - ptcl.pt_y_, scinti_front_intersec(2) - ptcl.pt_z_;
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
        std::cout << "void cilinder" << std::endl;
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
        std::cout << "no collision" << std::endl;
        return return_point;
    }
    ds = sqrt(ds);
    double a1 = (B - ds) / A,
        a2 = (B + ds) / A;
    return_point.at(2).at(0) = limtozero(ptcl.pt_x_ + a1 * v(0));
    return_point.at(2).at(1) = limtozero(ptcl.pt_y_ + a1 * v(1));
    return_point.at(2).at(2) = limtozero(ptcl.pt_z_ + a1 * v(2));
    return_point.at(3).at(0) = limtozero(ptcl.pt_x_ + a1 * v(0));
    return_point.at(3).at(1) = limtozero(ptcl.pt_y_ + a1 * v(1));
    return_point.at(3).at(2) = limtozero(ptcl.pt_z_ + a1 * v(2));

    std::cout << "side1: " << return_point.at(2).at(0) << ", " << return_point.at(2).at(1) << ", " << return_point.at(2).at(2) << std::endl;
    std::cout << "side2: " << return_point.at(3).at(0) << ", " << return_point.at(3).at(1) << ", " << return_point.at(3).at(2) << std::endl;
    
    return return_point;
}

double scinti::intersec_dist(particle initptcl, particle ptcl)
{
    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<double> zero_vector = { 0, 0, 0 };
    if (intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector && intersec.at(3) == zero_vector && !initpointcheck(initptcl, ptcl))
    {
        std::cout << "front only" << std::endl;
        return std::sqrt(std::pow(ptcl.pt_x_ - intersec.at(0).at(0), 2) + std::pow(ptcl.pt_y_ - intersec.at(0).at(1), 2) + std::pow(ptcl.pt_z_ - intersec.at(0).at(2), 2));
    }

    if (intersec.at(0) == zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector && intersec.at(3) == zero_vector && !initpointcheck(initptcl, ptcl))
    {
        std::cout << "back only" << std::endl;
        return std::sqrt(std::pow(ptcl.pt_x_-intersec.at(1).at(0), 2) + std::pow(ptcl.pt_y_ -intersec.at(1).at(1), 2) + std::pow(ptcl.pt_z_ -intersec.at(1).at(2), 2));
    }
    
    if (intersec.at(0) == zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector && intersec.at(3) == zero_vector && !initpointcheck(initptcl, ptcl))
    {
        std::cout << "side1 only" << std::endl;
        return std::sqrt(std::pow(ptcl.pt_x_-intersec.at(2).at(0), 2) + std::pow(ptcl.pt_y_ -intersec.at(2).at(1), 2) + std::pow(ptcl.pt_z_ -intersec.at(2).at(2), 2));
    }

    if (intersec.at(0) == zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector && intersec.at(3) != zero_vector && !initpointcheck(initptcl, ptcl))
    {
        std::cout << "side2 only" << std::endl;
        return std::sqrt(std::pow(ptcl.pt_x_-intersec.at(3).at(0), 2) + std::pow(ptcl.pt_y_ -intersec.at(3).at(1), 2) + std::pow(ptcl.pt_z_ -intersec.at(3).at(2), 2));
    }
    
    if (intersec.at(0) != zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector && intersec.at(3) == zero_vector)
    {
        std::cout << "front back" << std::endl;
        return std::sqrt(std::pow(intersec.at(0).at(0) - intersec.at(1).at(0), 2) + std::pow(intersec.at(0).at(1) - intersec.at(1).at(1), 2) + std::pow(intersec.at(0).at(2) - intersec.at(1).at(2), 2));
    }

    if (intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector && intersec.at(3) == zero_vector)
    {
        std::cout << "front side1" << std::endl;
        return std::sqrt(std::pow(intersec.at(0).at(0) - intersec.at(2).at(0), 2) + std::pow(intersec.at(0).at(1) - intersec.at(2).at(1), 2) + std::pow(intersec.at(0).at(2) - intersec.at(2).at(2), 2));
    }

    if (intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector && intersec.at(3) != zero_vector)
    {
        std::cout << "front side2" << std::endl;
        return std::sqrt(std::pow(intersec.at(0).at(0) - intersec.at(3).at(0), 2) + std::pow(intersec.at(0).at(1) - intersec.at(3).at(1), 2) + std::pow(intersec.at(0).at(2) - intersec.at(3).at(2), 2));
    }

    else
    {
        std::cout << "else" << std::endl;
        return 0;
    }
}