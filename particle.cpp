#include "particle.h"

double particle::limtozero(double num)
{
    if (0 < num && num < 0.00000000001)
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
    std::cout << "x: " << dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\ty: " << dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\tz: " << dist * std::cos(this->dir_theta_) << "\n" << "dist: " << sqrt(pow(dist * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), 2) + pow(dist * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), 2) + pow(dist * std::cos(this->dir_theta_), 2)) << std::endl;
}

// ‚Ü‚¾–¢Š®¬(turn_test‰ñ‚·‚Æ‚í‚©‚é)
void particle::turn(double angle)
{
    // std::cout << "loadeddeg: " << this->dir_theta_ << ", " << this->dir_phi_ << std::endl;
    std::random_device randseed_gen;
    std::mt19937 randengine(randseed_gen());
    std::uniform_real_distribution<> initrotate(0.0, 2.0);
    double randangle = initrotate(randengine)*M_PI;

    // ˆÈ‰º™R‚³‚ñ•û®
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
    // nt = t + b;
    // std::cout << "nt(0):" << nt(0) << " nt(1):" << nt(1) << " nt(2):" << nt(2) << std::endl;
    // this->dir_theta_ += atan(nt(1)/nt(0));
    // this->dir_phi_ += atan(sqrt(pow(nt(0), 2)+pow(nt(1), 2))/nt(2));

    // ˆÈ‰º‰ñ“]s—ñ•û®
    Eigen::Vector3d t, added_t, rotated_t;
    Eigen::Matrix3d rot3d;
    t << limtozero(std::sin(this->dir_theta_) * std::cos(this->dir_phi_)), limtozero(std::sin(this->dir_theta_) * std::sin(this->dir_phi_)), limtozero(std::cos(this->dir_theta_));
    // t.normalize();
    // std::cout << "t: " << t(0) << ", " << t(1) << ", " << t(2) << std::endl;
    added_t << limtozero(std::sin(this->dir_theta_+angle) * std::cos(this->dir_phi_)), limtozero(std::sin(this->dir_theta_+angle) * std::sin(this->dir_phi_)), limtozero(std::cos(this->dir_theta_+angle));
    // added_t.normalize();
    // std::cout << "added_t: " << added_t(0) << ", " << added_t(1) << ", " << added_t(2) << std::endl;
    rot3d = Eigen::AngleAxisd(randangle, t);
    rotated_t = rot3d * added_t;
    // rotated_t.normalize();
    // std::cout << "rotated_t: " << rotated_t(0) << ", " << rotated_t(1) << ", " << rotated_t(2) << std::endl;
    this->dir_theta_ += atan(rotated_t(1)/rotated_t(0));
    this->dir_phi_ += atan(sqrt(pow(rotated_t(0), 2)+pow(rotated_t(1), 2))/rotated_t(2));

    // ˆÈ‰º‚»‚Ì‚Ü‚Ü‘«‚·•û®
    // this->dir_theta_ += angle;
}

void particle::turn_test(double theta, double phi, double angle)
{
    std::ofstream turn_rect_angle("./data/turn_rect_angle.dat");
    turn_rect_angle << this->pt_x_ << "\t" << this->pt_y_ << "\t" << this->pt_z_ << "\n";
    this->dir_theta_ = theta;
    this->dir_phi_ = phi;
    turn_rect_angle << std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\t" << std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\t" << std::cos(this->dir_theta_) << "\n";
    
    for (int i = 0; i < 2500; i++)
    {
        this->dir_theta_ = theta;
        this->dir_phi_ = phi;
        turn(angle);
        turn_rect_angle << std::sin(this->dir_theta_) * std::cos(this->dir_phi_) << "\t" << std::sin(this->dir_theta_) * std::sin(this->dir_phi_) << "\t" << std::cos(this->dir_theta_) << "\n";
    }
    
}