#include "scinti.h"

#define MEC2 510.99895 // MeV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23) // mol^-1

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
    std::vector<std::vector<double>> return_point = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
    Eigen::Vector3d traject, scinti_centerline, scinti_frontcenter, scinti_front_intersec, scinti_backcenter, scinti_back_intersec;
    bool frontflag = false;

    traject << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
    traject.normalize();

    scinti_centerline << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_)* std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    scinti_centerline.normalize();
    std::cout << "scinti_ceneterline: " << scinti_centerline(0) << ", " << scinti_centerline(1) << ", " << scinti_centerline(2) << std::endl;

    scinti_frontcenter << this->pt_x_ + (this->depth_/2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (this->depth_ / 2) * std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (this->depth_ / 2) * std::cos(this->dir_theta_);
    std::cout << "scinti_frontcenter: " << scinti_frontcenter(0) << ", " << scinti_frontcenter(1) << ", " << scinti_frontcenter(2) << std::endl;

    double front_t = (scinti_centerline(0)*(scinti_frontcenter(0)-ptcl.pt_x_) + scinti_centerline(1) * (scinti_frontcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_frontcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (std::sqrt(std::pow(scinti_frontcenter(0) - (ptcl.pt_x_ + traject(0)*front_t), 2) + std::pow(scinti_frontcenter(1) - (ptcl.pt_y_ + traject(1) * front_t), 2) + std::pow(scinti_frontcenter(2) - (ptcl.pt_z_ + traject(2) * front_t), 2)) < this->rad_)
    {
        return_point.at(0).at(0) = ptcl.pt_x_ + traject(0) * front_t;
        return_point.at(0).at(1) = ptcl.pt_y_ + traject(1) * front_t;
        return_point.at(0).at(2) = ptcl.pt_z_ + traject(2) * front_t;
        frontflag = true;
        std::cout << "front: " << return_point.at(0).at(0) << ", " << return_point.at(0).at(1) << ", " << return_point.at(0).at(2) << std::endl;
    }

    scinti_backcenter << this->pt_x_ + (-this->depth_ / 2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (-this->depth_ / 2)* std::sin(this->dir_theta_)* std::sin(this->dir_phi_), this->pt_z_ + (-this->depth_ / 2)* std::cos(this->dir_theta_);
    std::cout << "scinti_backcenter: " << scinti_backcenter(0) << ", " << scinti_backcenter(1) << ", " << scinti_backcenter(2) << std::endl;
    double back_t = (scinti_centerline(0) * (scinti_backcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_backcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_backcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    if (std::sqrt(std::pow(scinti_backcenter(0) - (ptcl.pt_x_ + traject(0) * back_t), 2) + std::pow(scinti_backcenter(1) - (ptcl.pt_y_ + traject(1) * back_t), 2) + std::pow(scinti_backcenter(2) - (ptcl.pt_z_ + traject(2) * back_t), 2)) < this->rad_)
    {
        return_point.at(1).at(0) = ptcl.pt_x_ + traject(0) * back_t;
        return_point.at(1).at(1) = ptcl.pt_y_ + traject(1) * back_t;
        return_point.at(1).at(2) = ptcl.pt_z_ + traject(2) * back_t;
        std::cout << "back: " << return_point.at(1).at(0) << ", " << return_point.at(1).at(1) << ", " << return_point.at(1).at(2) << std::endl;
        return return_point;
    }

    else if(frontflag == true)
    {
        double move_t = front_t, dist = std::sqrt(2) * this->depth_;
        // std::cout << "move_t: " << move_t << std::endl;
        Eigen::Vector3d move_point, u, v, w;
        while (std::abs(dist - this->rad_) > 0.01)
        {
            // std::cout << "whileflag: " << std::abs(dist - this->rad_) << std::endl;
            if(move_t > 0)
            {
                move_t += 0.001;
            }
            else
            {
                move_t -= 0.001;
            }
            move_point << traject(0) * move_t, traject(1) * move_t, traject(2) * move_t;
            u << scinti_backcenter(0) - scinti_frontcenter(0), scinti_backcenter(1) - scinti_frontcenter(1), scinti_backcenter(2) - scinti_frontcenter(2);
            v << move_point(0) - scinti_frontcenter(0), move_point(1) - scinti_frontcenter(1), move_point(2) - scinti_frontcenter(2);
            w = u.cross(v)/u.norm();
            dist = w.norm();
            // std::cout << "dist: " << dist << std::endl;
            // std::cout << "move_point: " << move_point(0) << ", " << move_point(1) << ", " << move_point(2) << std::endl;
        }
        return_point.at(2).at(0) = move_point(0);
        return_point.at(2).at(1) = move_point(1);
        return_point.at(2).at(2) = move_point(2);
        std::cout << "side: " << return_point.at(2).at(0) << ", " << return_point.at(2).at(1) << ", " << return_point.at(2).at(2) << std::endl;
        return return_point;
    }
    return return_point;
}

double scinti::intersec_dist(particle ptcl)
{
    std::vector<std::vector<double>> intersec = scinti::intersec(ptcl);
    std::vector<double> zero_vector = { 0, 0, 0 };
    if ((intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) == zero_vector)||(intersec.at(0) == zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector)||(intersec.at(0) == zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector))
    {
        return std::max({std::sqrt(std::pow(ptcl.pt_x_ - intersec.at(0).at(0), 2) + std::pow(ptcl.pt_y_ - intersec.at(0).at(1), 2) + std::pow(ptcl.pt_z_ - intersec.at(0).at(2), 2)), std::sqrt(std::pow(ptcl.pt_x_-intersec.at(1).at(0), 2) + std::pow(ptcl.pt_y_ -intersec.at(1).at(1), 2) + std::pow(ptcl.pt_z_ -intersec.at(1).at(2), 2)), std::sqrt(std::pow(ptcl.pt_x_ - intersec.at(2).at(0), 2) + std::pow(ptcl.pt_y_ - intersec.at(2).at(1), 2) + std::pow(ptcl.pt_z_ - intersec.at(2).at(2), 2))});
    }

    else if (intersec.at(0) != zero_vector && intersec.at(1) != zero_vector && intersec.at(2) == zero_vector)
    {
        return std::sqrt(std::pow(intersec.at(0).at(0) - intersec.at(1).at(0), 2) + std::pow(intersec.at(0).at(1) - intersec.at(1).at(1), 2) + std::pow(intersec.at(0).at(2) - intersec.at(1).at(2), 2));
    }

    else if (intersec.at(0) != zero_vector && intersec.at(1) == zero_vector && intersec.at(2) != zero_vector)
    {
        return std::sqrt(std::pow(intersec.at(0).at(0) - intersec.at(2).at(0), 2) + std::pow(intersec.at(0).at(1) - intersec.at(2).at(1), 2) + std::pow(intersec.at(0).at(2) - intersec.at(2).at(2), 2));
    }
    else
    
    {
        return 0;
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