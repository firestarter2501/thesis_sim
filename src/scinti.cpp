#include "scinti.h"

#define MEC2 510.99895                         // keV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23)     // mol^-1

using vsize_t = std::vector<int>::size_type;

std::vector<bool> scinti::intersecenablecheck(std::vector<std::vector<double>> intersecvec, particle ptcl)
{
    std::vector<bool> returnvec;
    for (vsize_t i = 0; i < intersecvec.size(); i++)
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
        // std::cout << "dirnorm.norm(): " << dirnorm.norm() << std::endl;
        if (std::abs(centerdist.norm()) <= std::sqrt(2) * this->rad_ && dirnorm.norm() < 0.000001)
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

    if (std::abs(centerdist.norm()) <= std::sqrt(2) * this->rad_ && ec_truecount != 2)
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
        std::cout << "ptclinsidecheck: " << std::min({std::abs(traject1vec.norm()), std::abs(traject2vec.norm())}) << std::endl;
        return std::min({std::abs(traject1vec.norm()), std::abs(traject2vec.norm())});
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
    for (vsize_t i = 0; i < intersec.size(); i++)
    {
        std::cout << showfacetype(i) << ": ";
        for (vsize_t j = 0; j < intersec.at(0).size(); j++)
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

void scinti::initscinti(double pt_x, double pt_y, double pt_z, double theta, double phi, double depth, double z, double dens, double atomweight, double pmtsdevslope, double pmtsdevintersec)
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
    this->pmtsdevslope = pmtsdevslope;
    this->pmtsdevintersec = pmtsdevintersec;
    this->ene_buffer_ = 0;
}

std::vector<std::vector<double>> scinti::intersec(particle ptcl)
{
    std::vector<std::vector<double>> return_point = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Eigen::Vector3d traject, scinti_centerline, scinti_frontcenter, scinti_front_intersec, frontdist, scinti_backcenter, scinti_back_intersec, backdist;

    traject << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
    traject.normalize();

    // std::cout << "dir_theta_: " << this->dir_theta_ << ", dir_phi_: " << this->dir_phi_ << std::endl;

    scinti_centerline << std::sin(this->dir_theta_) * std::cos(this->dir_phi_), std::sin(this->dir_theta_) * std::sin(this->dir_phi_), std::cos(this->dir_theta_);
    scinti_centerline.normalize();
    // std::cout << "scinti_ceneterline: " << scinti_centerline(0) << ", " << scinti_centerline(1) << ", " << scinti_centerline(2) << std::endl;

    scinti_frontcenter << this->pt_x_ + (this->depth_ / 2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (this->depth_ / 2) * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), this->pt_z_ + (this->depth_ / 2) * std::cos(this->dir_theta_);
    // std::cout << "scinti_frontcenter: " << scinti_frontcenter(0) << ", " << scinti_frontcenter(1) << ", " << scinti_frontcenter(2) << std::endl;

    double front_t = (scinti_centerline(0) * (scinti_frontcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_frontcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_frontcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    scinti_front_intersec << ptcl.pt_x_ + traject(0) * front_t, ptcl.pt_y_ + traject(1) * front_t, ptcl.pt_z_ + traject(2) * front_t;
        
    return_point.at(0).at(0) = scinti_front_intersec(0);
    return_point.at(0).at(1) = scinti_front_intersec(1);
    return_point.at(0).at(2) = scinti_front_intersec(2);
    // std::cout << "front: " << return_point.at(0).at(0) << ", " << return_point.at(0).at(1) << ", " << return_point.at(0).at(2) << std::endl;

    scinti_backcenter << this->pt_x_ + (-this->depth_ / 2) * std::sin(this->dir_theta_) * std::cos(this->dir_phi_), this->pt_y_ + (-this->depth_ / 2) * std::sin(this->dir_theta_) * std::sin(this->dir_phi_), this->pt_z_ + (-this->depth_ / 2) * std::cos(this->dir_theta_);
    // std::cout << "scinti_backcenter: " << scinti_backcenter(0) << ", " << scinti_backcenter(1) << ", " << scinti_backcenter(2) << std::endl;

    double back_t = (scinti_centerline(0) * (scinti_backcenter(0) - ptcl.pt_x_) + scinti_centerline(1) * (scinti_backcenter(1) - ptcl.pt_y_) + scinti_centerline(2) * (scinti_backcenter(2) - ptcl.pt_z_)) / ((scinti_centerline(0) * traject(0)) + (scinti_centerline(1) * traject(1)) + (scinti_centerline(2) * traject(2)));
    scinti_back_intersec << ptcl.pt_x_ + traject(0) * back_t, ptcl.pt_y_ + traject(1) * back_t, ptcl.pt_z_ + traject(2) * back_t;

    return_point.at(1).at(0) = scinti_back_intersec(0);
    return_point.at(1).at(1) = scinti_back_intersec(1);
    return_point.at(1).at(2) = scinti_back_intersec(2);
    // std::cout << "back: " << return_point.at(1).at(0) << ", " << return_point.at(1).at(1) << ", " << return_point.at(1).at(2) << std::endl;

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
           rr = this->rad_ * this->rad_;
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

double scinti::scintillation(std::string scintidata, std::string csdata, particle &ptcl)
{
    std::ifstream initscinticonf("../data/" + scintidata + ".conf");
    std::vector<double> initscinticonf_list;
    std::string initline;
    while (std::getline(initscinticonf, initline))
    {
        initscinticonf_list.push_back(std::stod(initline));
        // std::cout << initline << std::endl;
    }
    initcs("../data/" + csdata + ".conf", this->crosssec_table_);
    this->initscinti(initscinticonf_list.at(0), initscinticonf_list.at(1), initscinticonf_list.at(2), initscinticonf_list.at(3) * M_PI, initscinticonf_list.at(4) * M_PI, initscinticonf_list.at(5), initscinticonf_list.at(6), initscinticonf_list.at(7), initscinticonf_list.at(8), initscinticonf_list.at(9), initscinticonf_list.at(10));

    bool react_flag = true;
    int react_count = 0;
    this->ene_buffer_ = 0;
    while (react_flag == true)
    {
        double traject_dist = this->intersec_dist(ptcl),
               pe_cs = crosssec(ptcl.ene_, 1, this->crosssec_table_),
               pe_len = reactlen(pe_cs, this->dens_),
               cs_ang = cs_angle(ptcl.ene_),
               cs_cs = crosssec(ptcl.ene_, 2, this->crosssec_table_),
               cs_len = reactlen(cs_cs, this->dens_),
               pp_cs = crosssec(ptcl.ene_, 3, this->crosssec_table_),
               pp_len = reactlen(pp_cs, this->dens_);

        if (0 < traject_dist && traject_dist < std::min({pe_len, cs_len, pp_len}) && ptcl.ene_ > 0)
        {
            react_flag = false;
            ptcl.move(this->ptclinsidecheck(ptcl) + traject_dist);
            std::cout << "too short(traject_dist: " << traject_dist << ", pe_len: " << pe_len << ", cs_len: " << cs_len << ", pp_len: " << pp_len << std::endl;
            std::cout << "ene_buffer_: " << this->ene_buffer_ << std::endl;
        }
        else if (traject_dist <= 0 || ptcl.ene_ <= 0)
        {
            react_flag = false;
            std::cout << "outside or zero ene(traject_dist: " << traject_dist << ", pe_len: " << pe_len << ", cs_len: " << cs_len << ", pp_len: " << pp_len << std::endl;
            std::cout << "ene_buffer_: " << this->ene_buffer_ << std::endl;
        }
        else
        {
            react_count++;
            std::cout << "before ene_buffer_: " << this->ene_buffer_ << std::endl;
            if (pe_len <= cs_len && pe_len <= pp_len)
            {
                this->ene_buffer_ += normdist(ptcl.ene_, lineareq(ptcl.ene_, this->pmtsdevslope, this->pmtsdevintersec));
                ptcl.ene_ = 0;
                react_flag = false;
                std::cout << "---pe---" << std::endl;
                std::cout << "ene_buffer: " << this->ene_buffer_ << std::endl;
            }
            else if (cs_len <= pe_len && cs_len <= pp_len)
            {
                // Eigen::Vector3d beforepoint, afterpoint, moveddist;
                this->ene_buffer_ += normdist(ptcl.ene_ - scatphotonene(ptcl.ene_, cs_ang), lineareq(ptcl.ene_ - scatphotonene(ptcl.ene_, cs_ang), this->pmtsdevslope, this->pmtsdevintersec));
                ptcl.ene_ = scatphotonene(ptcl.ene_, cs_ang);
                // beforepoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                ptcl.move(this->ptclinsidecheck(ptcl) + cs_len);
                // afterpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                // moveddist = afterpoint - beforepoint;
                // std::cout << "moved dist: " << std::abs(moveddist.norm()) << std::endl;
                ptcl.turn(cs_ang);
                react_flag = true;
                std::cout << "---cs---" << std::endl;
                std::cout << "ene_buffer_: " << this->ene_buffer_ << " cs_len: " << cs_len << " cs_ang: " << cs_ang << std::endl;
            }
            else if (pp_len <= pe_len && pp_len <= cs_len)
            {
                // Eigen::Vector3d beforepoint, afterpoint, moveddist;
                ptcl.ene_ = MEC2;
                ptcl.initptcl(ptcl.ene_, ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_);
                // beforepoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                ptcl.move(this->ptclinsidecheck(ptcl) + pp_len);
                // afterpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                // moveddist = afterpoint - beforepoint;
                // std::cout << "moved dist: " << std::abs(moveddist.norm()) << std::endl;
                react_flag = true;
                std::cout << "---pp---" << std::endl;
                std::cout << "ene_buffer_: " << this->ene_buffer_ << " pp_len: " << pp_len << std::endl;
            }
            else
            {
                std::cout << "judge error" << std::endl;
                react_flag = false;
            }
        }

        if (!react_flag && 0 < this->ene_buffer_)
        {
                return this->ene_buffer_;
                std::cout << "sum_ene check&add done" << std::endl;
        }
    }
    // if (0 < ptcl.ene_ && !absorp_flag)
    // {
    //     react_flag = true;
    // }
    return -1;
}