#pragma once
#include <random>
#include <fstream>
#include "scinti.h"

void crosssec_test(scinti scintillator);

// Œõ“d‹zû’f–ÊÏ‚ğ•Ô‚·ŠÖ”
double pe_crosssec(double ene, int z);

// ƒRƒ“ƒvƒgƒ“U—‚·‚é‚Æ‚«‚ÌU—Œõq‚ÌƒGƒlƒ‹ƒM[‚ğ•Ô‚·ŠÖ”
double scatphotonene(double ene, double angle);

// ƒNƒ‰ƒCƒ“m‰È‚Ì®‚ğŒvZ‚·‚éŠÖ”
double kleinnishinaeq(double ene, double angle);

// von Neumann‚ÌŠü‹p–@‚ÅƒRƒ“ƒvƒgƒ“U—‚ÌŠp“x•ª•z‚ÉŠî‚Ã‚¢‚½U—Šp‚Ì—”‚ğ•Ô‚·ŠÖ”
double cs_angle(double ene);

void cs_angle_test(int num);

double reactlen(double crosssec, double dens);

double normdist(double mean, double sdev);

int randchoice(std::vector<double> probvec);

double lineareq(double ene, double slope, double intersec);

void showinfo(std::vector<particle> photon, double traject_dist, double pe_len, double cs_len, double pp_len);