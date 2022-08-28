#include "scinti.h"

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * pow(10, -13) // cm
#define NUMA 6.02214076 * pow(10, 23) // mol^-1

void scinti::initsinti(double x, double y, double z, double theta, double phi, double depth, double dens, double atomweight)
{
    this->pt_x_ = x;
    this->pt_y_ = y;
    this->pt_z_ = z;
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
    std::vector<std::vector<double>> return_point = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
    Eigen::Vector3d traject, scinti_centerline, scinti_frontcenter, scinti_front_intersec, scinti_backcenter, scinti_back_intersec;

    traject << sin(ptcl.dir_theta_) * cos(ptcl.dir_phi_), sin(ptcl.dir_theta_) * sin(ptcl.dir_phi_), cos(ptcl.dir_theta_);
    traject.normalize();

    scinti_centerline << sin(this->dir_theta_) * cos(this->dir_phi_), sin(this->dir_theta_)* sin(this->dir_phi_), cos(this->dir_theta_);
    scinti_centerline.normalize();

    scinti_frontcenter << (this->depth_/2) * sin(this->dir_theta_) * cos(this->dir_phi_), (this->depth_ / 2) * sin(this->dir_theta_)* sin(this->dir_phi_), (this->depth_ / 2) * cos(this->dir_theta_);
    double front_t = (scinti_centerline(0)*(scinti_frontcenter(0)-ptcl.pt_x_) + scinti_centerline(1) * (scinti_frontcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_frontcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (sqrt(pow(scinti_frontcenter(0) - (ptcl.pt_x_ + traject(0)*front_t), 2) + pow(scinti_frontcenter(1) - (ptcl.pt_y_ + traject(1) * front_t), 2) + pow(scinti_frontcenter(2) - (ptcl.pt_z_ + traject(2) * front_t), 2)) < this->rad_)
    {
        return_point.at(0).at(0) = ptcl.pt_x_ + traject(0) * front_t;
        return_point.at(0).at(1) = ptcl.pt_y_ + traject(1) * front_t;
        return_point.at(0).at(2) = ptcl.pt_z_ + traject(2) * front_t;
    }

    scinti_backcenter << (-this->depth_ / 2) * sin(this->dir_theta_) * cos(this->dir_phi_), (-this->depth_ / 2)* sin(this->dir_theta_)* sin(this->dir_phi_), (-this->depth_ / 2)* cos(this->dir_theta_);
    double back_t = (scinti_centerline(0) * (scinti_backcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_backcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_backcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (sqrt(pow(scinti_backcenter(0) - (ptcl.pt_x_ + traject(0) * back_t), 2) + pow(scinti_backcenter(1) - (ptcl.pt_y_ + traject(1) * back_t), 2) + pow(scinti_backcenter(2) - (ptcl.pt_z_ + traject(2) * back_t), 2)) < this->rad_)
    {
        return_point.at(1).at(0) = ptcl.pt_x_ + traject(0) * back_t;
        return_point.at(1).at(1) = ptcl.pt_y_ + traject(1) * back_t;
        return_point.at(1).at(2) = ptcl.pt_z_ + traject(2) * back_t;
        
        return return_point;
    }

    double start_t = front_t, end_t = front_t + sqrt(2) * this->depth_, dist = sqrt(2) * this->depth_;
    Eigen::Vector3d move_point, u, v, w;
    while (sqrt(pow(dist - this->rad_, 2)) < 0.01)
    {
        move_point << ptcl.pt_x_ + traject(0) * (end_t - start_t)/2, ptcl.pt_y_ + traject(1) * (end_t - start_t)/2, ptcl.pt_z_ + traject(2) * (end_t - start_t)/2;
        u << scinti_backcenter(0) - scinti_frontcenter(0), scinti_backcenter(1) - scinti_frontcenter(1), scinti_backcenter(2) - scinti_frontcenter(2);
        v << move_point(0) - scinti_frontcenter(0), move_point(1) - scinti_frontcenter(1), move_point(2) - scinti_frontcenter(2);
        w = u.cross(v)/u.norm();
        dist = w.norm();
        if (dist < this->rad_)
        {
            end_t = (end_t - start_t) / 2;
        }
        else
        {
            start_t = (end_t - start_t) / 2;
        }
    }
    return_point.at(2).at(0) = move_point(0);
    return_point.at(2).at(1) = move_point(1);
    return_point.at(2).at(2) = move_point(2);
    return return_point;
}

double scinti::intersec_dist(particle ptcl)
{
    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<double> zero_vector = { 0, 0, 0 };
    if ((intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector)||(intersec.at(0) == zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector)||(intersec.at(0) == zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector))
    {
        return std::max(sqrt(pow(ptcl.pt_x_ - intersec.at(0).at(0), 2) + pow(ptcl.pt_y_ - intersec.at(0).at(1), 2) + pow(ptcl.pt_z_ - intersec.at(0).at(2), 2)), sqrt(pow(ptcl.pt_x_-intersec.at(1).at(0), 2) + pow(ptcl.pt_y_ -intersec.at(1).at(1), 2) + pow(ptcl.pt_z_ -intersec.at(1).at(2), 2)), sqrt(pow(ptcl.pt_x_ - intersec.at(2).at(0), 2) + pow(ptcl.pt_y_ - intersec.at(2).at(1), 2) + pow(ptcl.pt_z_ - intersec.at(2).at(2), 2)));
    }
    else
    {
        return std::max(sqrt(pow(intersec.at(0).at(0) - intersec.at(1).at(0), 2) + pow(intersec.at(0).at(1) - intersec.at(1).at(1), 2) + pow(intersec.at(0).at(2) - intersec.at(1).at(2), 2)), sqrt(pow(intersec.at(0).at(0) - intersec.at(2).at(0), 2) + pow(intersec.at(0).at(1) - intersec.at(2).at(1), 2) + pow(intersec.at(0).at(2) - intersec.at(2).at(2), 2)));
    }
}

bool scinti::internal_judge(particle ptcl)
{
    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<double> zero_vector = { 0, 0, 0 };
    if ((intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector) || (intersec.at(0) == zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector) || (intersec.at(0) == zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector))
    {
        return true;
    }
    else
    {
        return false;
    }
}