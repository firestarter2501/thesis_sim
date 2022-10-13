#include "block.h"

#define MEC2 510.99895 // keV
#define RELEC 2.8179403262 * std::pow(10, -13) // cm
#define NUMA 6.02214076 * std::pow(10, 23) // mol^-1

using vsize_t = std::vector<int>::size_type;

void block::initblock(double pt_x, double pt_y, double pt_z, double height, double width, double depth, double z, double dens, double atomweight)
{
    this->pt_x_ = pt_x;
    this->pt_y_ = pt_y;
    this->pt_z_ = pt_z;
    this->height_ = height;
    this->width_ = width;
    this->depth_ = depth;
    this->z_ = z;
    this->dens_ = dens;
    this->ndens_ = (dens / atomweight) * NUMA;
    this->atomweight_ = atomweight;
}

std::vector<std::vector<double>> block::intersec(particle ptcl)
{
    std::vector<std::vector<double>> return_point;

    Eigen::Vector3d traject;

    traject << std::sin(ptcl.dir_theta_) * std::cos(ptcl.dir_phi_), std::sin(ptcl.dir_theta_) * std::sin(ptcl.dir_phi_), std::cos(ptcl.dir_theta_);
    traject.normalize();

    Eigen::Vector3d x_centerline, front_center, front_intersec, back_center, back_intersec;
    x_centerline << 1, 0, 0;
    front_center << this->pt_x_+(this->depth_/2), this->pt_y_, this->pt_z_;
    double front_t = (x_centerline(0) * (front_center(0) - ptcl.pt_x_) + x_centerline(1) * (front_center(1) - ptcl.pt_y_) + x_centerline(2) * (front_center(2) - ptcl.pt_z_)) / ((x_centerline(0) * traject(0)) + (x_centerline(1) * traject(1)) + (x_centerline(2) * traject(2)));
    front_intersec << ptcl.pt_x_ + traject(0) * front_t, ptcl.pt_y_ + traject(1) * front_t, ptcl.pt_z_ + traject(2) * front_t;
    return_point.push_back({front_intersec(0), front_intersec(1), front_intersec(2)});
    
    back_center << this->pt_x_-(this->depth_/2), this->pt_y_, this->pt_z_;
    double back_t = (x_centerline(0) * (back_center(0) - ptcl.pt_x_) + x_centerline(1) * (back_center(1) - ptcl.pt_y_) + x_centerline(2) * (front_center(2) - ptcl.pt_z_)) / ((x_centerline(0) * traject(0)) + (x_centerline(1) * traject(1)) + (x_centerline(2) * traject(2)));
    back_intersec << ptcl.pt_x_ + traject(0) * back_t, ptcl.pt_y_ + traject(1) * back_t, ptcl.pt_z_ + traject(2) * back_t;
    return_point.push_back({back_intersec(0), back_intersec(1), back_intersec(2)});

    Eigen::Vector3d y_centerline, left_center, left_intersec, right_center, right_intersec;
    y_centerline << 0, 1, 0;
    left_center << this->pt_x_, this->pt_y_-(this->width_/2), this->pt_z_;
    double left_t = (y_centerline(0) * (left_center(0) - ptcl.pt_x_) + y_centerline(1) * (left_center(1) - ptcl.pt_y_) + y_centerline(2) * (left_center(2) - ptcl.pt_z_)) / ((y_centerline(0) * traject(0)) + (y_centerline(1) * traject(1)) + (y_centerline(2) * traject(2)));
    left_intersec << ptcl.pt_x_ + traject(0) * left_t, ptcl.pt_y_ + traject(1) * left_t, ptcl.pt_z_ + traject(2) * left_t;
    return_point.push_back({left_intersec(0), left_intersec(1), left_intersec(2)});
    
    right_center << this->pt_x_, this->pt_y_+(this->width_/2), this->pt_z_;
    double right_t = (y_centerline(0) * (right_center(0) - ptcl.pt_x_) + y_centerline(1) * (right_center(1) - ptcl.pt_y_) + y_centerline(2) * (right_center(2) - ptcl.pt_z_)) / ((y_centerline(0) * traject(0)) + (y_centerline(1) * traject(1)) + (y_centerline(2) * traject(2)));
    right_intersec << ptcl.pt_x_ + traject(0) * right_t, ptcl.pt_y_ + traject(1) * right_t, ptcl.pt_z_ + traject(2) * right_t;
    return_point.push_back({right_intersec(0), right_intersec(1), right_intersec(2)});

    Eigen::Vector3d z_centerline, top_center, top_intersec, bottom_center, bottom_intersec;
    z_centerline << 0, 0, 1;
    top_center << this->pt_x_, this->pt_y_, this->pt_z_+(this->height_/2);
    double top_t = (z_centerline(0) * (top_center(0) - ptcl.pt_x_) + z_centerline(1) * (top_center(1) - ptcl.pt_y_) + z_centerline(2) * (top_center(2) - ptcl.pt_z_)) / ((z_centerline(0) * traject(0)) + (z_centerline(1) * traject(1)) + (z_centerline(2) * traject(2)));
    top_intersec << ptcl.pt_x_ + traject(0) * top_t, ptcl.pt_y_ + traject(1) * top_t, ptcl.pt_z_ + traject(2) * top_t;
    return_point.push_back({top_intersec(0), top_intersec(1), top_intersec(2)});
    
    bottom_center << this->pt_x_, this->pt_y_, this->pt_z_-(this->height_/2);
    double bottom_t = (z_centerline(0) * (bottom_center(0) - ptcl.pt_x_) + z_centerline(1) * (bottom_center(1) - ptcl.pt_y_) + z_centerline(2) * (bottom_center(2) - ptcl.pt_z_)) / ((z_centerline(0) * traject(0)) + (z_centerline(1) * traject(1)) + (z_centerline(2) * traject(2)));
    bottom_intersec << ptcl.pt_x_ + traject(0) * bottom_t, ptcl.pt_y_ + traject(1) * bottom_t, ptcl.pt_z_ + traject(2) * bottom_t;
    return_point.push_back({bottom_intersec(0), bottom_intersec(1), bottom_intersec(2)});

    return return_point;
}

std::vector<bool> block::intersecenablecheck(std::vector<std::vector<double>> intersecvec, particle ptcl)
{
    std::vector<bool> returnvec;
    for (vsize_t i = 0; i < intersecvec.size(); i++)
    {
        if ((this->pt_x_-(this->depth_/2) <= intersecvec.at(i).at(0) && intersecvec.at(i).at(0) <= this->pt_x_+(this->depth_/2)) && (this->pt_y_-(this->width_/2) <= intersecvec.at(i).at(1) && intersecvec.at(i).at(1) <= this->pt_y_+(this->width_/2)) && (this->pt_z_-(this->height_/2) <= intersecvec.at(i).at(2) && intersecvec.at(i).at(2) <= this->pt_z_+(this->height_/2)))
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

double block::ptclinsidecheck(particle ptcl)
{
    Eigen::Vector3d ptclpoint;
    ptclpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
    std::vector<std::vector<double>> intersec = block::intersec(ptcl);
    std::vector<bool> enablecheck = intersecenablecheck(intersec, ptcl);
    int ec_truecount = std::count(enablecheck.begin(), enablecheck.end(), true);

    if (ec_truecount != 2)
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

std::string block::showfacetype(int type)
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
        return "left";
    }
    else if (type == 3)
    {
        return "right";
    }
    else if (type == 4)
    {
        return "top";
    }
    else if (type == 5)
    {
        return "bottom";
    }
    else
    {
        return "facetypeerror";
    }
}

double block::intersec_dist(particle ptcl)
{
    std::vector<std::vector<double>> intersec = block::intersec(ptcl);
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

void block::react(std::string blockdata, std::string csdata, particle &ptcl, bool &react_flag, bool &absorp_flag)
{
    std::ifstream blockdataconf("../data/" + blockdata + ".conf");
    std::vector<double> blockdataconf_list;
    std::string initline;
    while (std::getline(blockdataconf, initline))
    {
        blockdataconf_list.push_back(std::stod(initline));
        // std::cout << initline << std::endl;
    }
    initcs("../data/" + csdata + ".conf", this->crosssec_table_);
    this->initblock(blockdataconf_list.at(0), blockdataconf_list.at(1), blockdataconf_list.at(2), blockdataconf_list.at(3), blockdataconf_list.at(4), blockdataconf_list.at(5), blockdataconf_list.at(6), blockdataconf_list.at(7), blockdataconf_list.at(8));

    react_flag = true;
    int react_count = 0;
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

        if (traject_dist < std::min({pe_len, cs_len, pp_len}) || traject_dist == -1 /* || (scinti_num == 1 && photon.back().ene_ == ray_list.back().ene_)*/ /* || react_count >= 1*/)
        {
            react_flag = false;
            std::cout << "outside or too short(traject_dist: " << traject_dist << ", pe_len: " << pe_len << ", cs_len: " << cs_len << ", pp_len: " << pp_len << std::endl;
        }
        else
        {
            react_count++;
            if (pe_len <= cs_len && pe_len <= pp_len)
            {
                react_flag = false;
                absorp_flag = true;
                std::cout << "---pe---" << std::endl;
            }
            else if (cs_len <= pe_len && cs_len <= pp_len)
            {
                // Eigen::Vector3d beforepoint, afterpoint, moveddist;
                ptcl.ene_ = scatphotonene(ptcl.ene_, cs_ang);
                // beforepoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                ptcl.move(this->ptclinsidecheck(ptcl) + cs_len);
                // afterpoint << ptcl.pt_x_, ptcl.pt_y_, ptcl.pt_z_;
                // moveddist = afterpoint - beforepoint;
                // std::cout << "moved dist: " << std::abs(moveddist.norm()) << std::endl;
                ptcl.turn(cs_ang);
                react_flag = true;
                std::cout << "---cs---" << std::endl;
                std::cout << " cs_len: " << cs_len << " cs_ang: " << cs_ang << std::endl;
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
                std::cout << " pp_len: " << pp_len << std::endl;
            }
            else
            {
                std::cout << "judge error" << std::endl;
                react_flag = false;
            }
        }
    }
    // if (0 < ptcl.ene_ && !absorp_flag)
    // {
    //     react_flag = true;
    // }
}