#include "particle.h"

void particle::initptcl(double ene, double x, double y, double z)
{
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> inittheta(-1.0, 1.0);
    std::uniform_real_distribution<> initphi(0.0, 2.0);
    this->ene_ = ene;
    this->pt_x_ = x;
    this->pt_y_ = y;
    this->pt_z_ = z;
    double theta = std::acos(inittheta(randengine));
    while (theta < 0 || M_PI < theta)
    {
        theta = std::acos(inittheta(randengine));
    }
    
    this->dir_theta_ = theta;
    this->dir_phi_ = initphi(randengine) * M_PI;
};

void particle::initptcl_test(int num)
{
    std::ofstream initptcl_polar_angle("initptcl_polar_angle.dat"), initptcl_rectangular_angle("initptcl_rectangular_angle.dat");

    for (int i = 0; i < num; i++)
    {
        this->initptcl(0, 0, 0, 0);
        initptcl_polar_angle << this->dir_phi_ << "\t" << this->dir_theta_ << "\n";
        initptcl_rectangular_angle << std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\t" << std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\t" << std::cos(this->dir_theta_) << "\n";
    }
}

void particle::move(double dist)
{
    this->pt_x_ += dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_);
    this->pt_y_ += dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_);
    this->pt_z_ += dist * std::cos(this->dir_theta_);
}

void particle::move_test(double dist)
{
    this->initptcl(0, 0, 0, 0);
    this->move(dist);
    std::cout << "x: " << dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\ty: " << dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\tz: " << dist * std::cos(this->dir_theta_) << "\n" << "dist: " << sqrt(pow(dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), 2) + pow(dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), 2) + pow(dist * std::cos(this->dir_theta_), 2)) << std::endl;
}

// ‚Ü‚¾–¢Š®¬(turn_test‰ñ‚·‚Æ‚í‚©‚é)
void particle::turn(double angle)
{
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> initrotate(0.0, 2.0);
    Eigen::Vector3d point_in, point_out;
    Eigen::Matrix3d AxisAngle;
    Eigen::Vector3d axis;
    double rotateangle = initrotate(randengine) * M_PI;
    double const d = 1000;

    point_in << d * std::sin(this->dir_theta_ - angle) * std::cos(this->dir_phi_), d* std::sin(this->dir_theta_ - angle)* std::sin(this->dir_phi_), d* std::cos(this->dir_theta_ - angle);
    axis << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_)* std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    axis.normalize();
    AxisAngle = Eigen::AngleAxisd(rotateangle, axis);
    point_out = AxisAngle * point_in;
    this->dir_theta_ += std::atan(point_out(1, 0) / point_out(0, 0));
    this->dir_phi_ += std::atan(sqrt(pow(point_out(0, 0), 2) + pow(point_out(1, 0), 2)) / point_out(2, 0));
}

void particle::turn_test(double angle)
{
    this->dir_theta_ = 0;
    this->dir_phi_ = 0;
    double dist = 1000;
    double before_theta = this->dir_theta_, before_phi = this->dir_phi_, before_x = dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), before_y = dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), before_z = dist * std::cos(this->dir_theta_);
    std::cout << "theta: " << before_theta << "\tphi: " << before_phi << "\ninput_angle: " << angle << std::endl;
    this->turn(angle);
    double after_x = dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), after_y = dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), after_z = dist * std::cos(this->dir_theta_);
    std::cout << "theta: " << this->dir_theta_ << "\tphi: " << this->dir_phi_ << "\nangle: " << std::acos(((before_x * after_x) + (before_y * after_y) + (before_z * after_z)) / (sqrt(pow(before_x, 2) + pow(before_y, 2) + pow(before_z, 2)) * sqrt(pow(after_x, 2) + pow(after_y, 2) + pow(after_z, 2)))) << std::endl;
}