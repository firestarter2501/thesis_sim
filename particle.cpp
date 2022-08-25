#include "particle.h"

using namespace std;
using namespace Eigen;

void particle::initptcl(double ene, double x, double y, double z)
{
    random_device randseed_gen;
    mt19937 randengine(randseed_gen());
    uniform_real_distribution<> inittheta(-1.0, 1.0);
    uniform_real_distribution<> initphi(0.0, 2.0);
    this->ene_ = ene;
    this->pt_x_ = x;
    this->pt_y_ = y;
    this->pt_z_ = z;
    this->dir_theta_ = acos(inittheta(randengine));
    this->dir_phi_ = initphi(randengine) * M_PI;
};

void particle::initptcl_test(int num)
{
    ofstream initptcl_polar_angle("initptcl_polar_angle.dat"), initptcl_rectangular_angle("initptcl_rectangular_angle.dat");

    for (int i = 0; i < num; i++)
    {
        this->initptcl(0, 0, 0, 0);
        initptcl_polar_angle << this->dir_phi_ << "\t" << this->dir_theta_ << "\n";
        initptcl_rectangular_angle << sin(this->dir_theta_) * cos(this->dir_phi_) << "\t" << sin(this->dir_theta_) * sin(this->dir_phi_) << "\t" << cos(this->dir_theta_) << "\n";
    }
}

void particle::move(double dist)
{
    this->pt_x_ += dist * sin(this->dir_theta_) * cos(this->dir_phi_);
    this->pt_y_ += dist * sin(this->dir_theta_) * sin(this->dir_phi_);
    this->pt_z_ += dist * cos(this->dir_theta_);
}

void particle::move_test(double dist)
{
    this->initptcl(0, 0, 0, 0);
    this->move(dist);
    cout << "x: " << dist * sin(this->dir_theta_) * cos(this->dir_phi_) << "\ty: " << dist * sin(this->dir_theta_) * sin(this->dir_phi_) << "\tz: " << dist * cos(this->dir_theta_) << "\n" << "dist: " << sqrt(pow(dist * sin(this->dir_theta_) * cos(this->dir_phi_), 2) + pow(dist * sin(this->dir_theta_) * sin(this->dir_phi_), 2) + pow(dist * cos(this->dir_theta_), 2)) << endl;
}

void particle::turn(double angle)
{
    random_device randseed_gen;
    mt19937 randengine(randseed_gen());
    uniform_real_distribution<> initrotate(0.0, 2.0);
    Vector3d point_in, point_out;
    Matrix3d AxisAngle;
    Vector3d axis;
    double rotateangle = initrotate(randengine) * M_PI;
    double const d = 1000;

    point_in << d * sin(this->dir_theta_ - angle) * cos(this->dir_phi_), d* sin(this->dir_theta_ - angle)* sin(this->dir_phi_), d* cos(this->dir_theta_ - angle);
    axis << sin(this->dir_theta_) * cos(this->dir_phi_), sin(this->dir_theta_)* sin(this->dir_phi_), cos(this->dir_theta_);
    axis.normalize();
    AxisAngle = AngleAxisd(rotateangle, axis);
    point_out = AxisAngle * point_in;
    this->dir_theta_ += atan(point_out(1, 0) / point_out(0, 0));
    this->dir_phi_ += atan(sqrt(pow(point_out(0, 0), 2) + pow(point_out(1, 0), 2)) / point_out(2, 0));
}

void particle::turn_test(double angle)
{
    this->dir_theta_ = 0;
    this->dir_phi_ = 0;
    double dist = 1000;
    double before_theta = this->dir_theta_, before_phi = this->dir_phi_, before_x = dist * sin(this->dir_theta_) * cos(this->dir_phi_), before_y = dist * sin(this->dir_theta_) * sin(this->dir_phi_), before_z = dist * cos(this->dir_theta_);
    cout << "theta: " << before_theta << "\tphi: " << before_phi << "\ninput_angle: " << angle << endl;
    this->turn(angle);
    double after_x = dist * sin(this->dir_theta_) * cos(this->dir_phi_), after_y = dist * sin(this->dir_theta_) * sin(this->dir_phi_), after_z = dist * cos(this->dir_theta_);
    cout << "theta: " << this->dir_theta_ << "\tphi: " << this->dir_phi_ << "\nangle: " << acos(((before_x * after_x) + (before_y * after_y) + (before_z * after_z)) / (sqrt(pow(before_x, 2) + pow(before_y, 2) + pow(before_z, 2)) * sqrt(pow(after_x, 2) + pow(after_y, 2) + pow(after_z, 2)))) << endl;
}