#include "particle.h"

double particle::limtozero(double num)
{
    if (-0.00000001 < num && num < 0.00000001)
    {
        return 0;
    }
    else
    {
        return num;
    }
}

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
//std::cout << "x: " << dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\ty: " << dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\tz: " << dist * std::cos(this->dir_theta_) << "\n" << "dist: " << sqrt(pow(dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), 2) + pow(dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), 2) + pow(dist * std::cos(this->dir_theta_), 2)) << std::endl;
}

// まだ未完成(turn_test回すとわかる)
void particle::turn(double angle)
{
    // std::cout << "loadeddeg: " << this->dir_theta_ << ", " << this->dir_phi_ << std::endl;
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> initrotate(0, 2);

    // 以下杉山さん方式
    // Eigen::Vector3d t, n, q, b, nt;
    // t << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_) * std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    // t.normalize();
    // std::cout << "t(0):" << t(0) << " t(1):" << t(1) << " t(2):" << t(2) << std::endl;
    // n << ((t(0)*t(1))/(pow(t(0), 2) + pow(t(1), 2))), ((t(1)*t(2))/(pow(t(0), 2) + pow(t(1), 2))), -sqrt(pow(t(0), 2) + pow(t(1), 2));
    // std::cout << "n(0):" << n(0) << " n(1):" << n(1) << " n(2):" << n(2) << std::endl;
    // q = n.cross(t);
    // std::cout << "q(0):" << q(0) << " q(1):" << q(1) << " q(2):" << q(2) << std::endl;
    // b << ((t(0)*t(2)*cos(randangle))/sqrt(pow(t(0), 2)+pow(t(1), 2)))-(t(1)*sin(randangle)), ((t(1)*t(2)*cos(randangle))/sqrt(pow(t(0), 2)+pow(t(1), 2)))+(t(1)*sin(randangle)), -sqrt(pow(t(0), 2)+pow(t(1), 2))*cos(randangle);
    // b.normalize();
    // std::cout << "b(0):" << b(0) << " b(1):" << b(1) << " b(2):" << b(2) << std::endl;
    // nt = angle*t + (M_PI-angle)*b;
    // std::cout << "nt(0):" << nt(0) << " nt(1):" << nt(1) << " nt(2):" << nt(2) << std::endl;
    // this->dir_theta_ += acos(nt(2)/nt.norm());
    // this->dir_phi_ += atan(nt(1)/nt(0));

    // 以下回転行列方式
    Eigen::Vector3d t, added_t, rotated_t;
    Eigen::Matrix3d rot3d;
    double distcheck, randangle;
    t << limtozero(std::sin(this->dir_theta_) * std::cos(this->dir_phi_)), limtozero(std::sin(this->dir_theta_) * std::sin(this->dir_phi_)), limtozero(std::cos(this->dir_theta_));
    t.normalize();
    added_t << limtozero(std::sin(this->dir_theta_-angle) * std::cos(this->dir_phi_)), limtozero(std::sin(this->dir_theta_-angle) * std::sin(this->dir_phi_)), limtozero(std::cos(this->dir_theta_-angle));
    randangle = initrotate(randengine)*M_PI;
    rotated_t = added_t*cos(randangle) + (1-cos(randangle))*(added_t.dot(t))*t + t.cross(added_t)*sin(randangle);
    rotated_t.normalize();
    distcheck = abs(t(0)*rotated_t(0)+t(1)*rotated_t(1)+t(2)*rotated_t(2))/t.norm();
    // std::cout << distcheck << "/" << cos(angle) << std::endl;
    if (limtozero(distcheck-cos(angle)) != 0)
    {
        Eigen::Vector3d t_flat, t_flat_mirror, zvec, t_vert1, t_vert2, t_vert1_mirror, t_vert2_mirror;
        t_flat = -t;
        t_flat(2) = 0;
        t_flat_mirror = rotated_t - 2*(rotated_t.dot(t_flat))*t_flat;
        zvec << 0, 0, 1;
        t_vert1 = t_flat*cos(this->dir_phi_+M_PI/4) + (1-cos(this->dir_phi_+M_PI/4))*(t_flat.dot(zvec))*zvec + zvec.cross(t_flat)*sin(this->dir_phi_+M_PI/4);
        t_vert1_mirror = t_flat_mirror - 2*(t_vert1.dot(t_flat_mirror))*t_flat_mirror;
        t_vert2 = t_flat*cos(this->dir_phi_-M_PI/4) + (1-cos(this->dir_phi_-M_PI/4))*(t_flat.dot(zvec))*zvec + zvec.cross(t_flat)*sin(this->dir_phi_-M_PI/4);
        t_vert2_mirror = t_flat_mirror - 2*(t_vert2.dot(t_flat_mirror))*t_flat_mirror;
        if (t_vert1_mirror(0) <= 0)
        {
            this->dir_theta_ = acos(t_vert1_mirror(2)/t_vert1_mirror.norm());
            this->dir_phi_ = atan(t_vert1_mirror(1)/t_vert1_mirror(0));
        }
        else
        {
            this->dir_theta_ = acos(t_vert2_mirror(2)/t_vert2_mirror.norm());
            this->dir_phi_ = atan(t_vert2_mirror(1)/t_vert2_mirror(0));
        }
    }
    else
    {
        this->dir_theta_ = acos(rotated_t(2)/rotated_t.norm());
        this->dir_phi_ = atan(rotated_t(1)/rotated_t(0));
    }
}

void particle::turn_test(double theta, double phi, double angle)
{
    std::ofstream turn_rect_angle("./data/turn_rect_angle.dat");
    turn_rect_angle << this->pt_x_ << "\t" << this->pt_y_ << "\t" << this->pt_z_ << "\n";
    this->dir_theta_ = theta;
    this->dir_phi_ = phi;
    // std::cout << "loadeddeg: " << this->dir_theta_ << ", " << this->dir_phi_ << std::endl;
    turn_rect_angle << std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\t" << std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\t" << std::cos(this->dir_theta_) << "\n";
    
    for (int i = 0; i < 2500; i++)
    {
        this->dir_theta_ = theta;
        this->dir_phi_ = phi;
        turn(angle);
        turn_rect_angle << std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\t" << std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\t" << std::cos(this->dir_theta_) << "\n";
    }
    
}