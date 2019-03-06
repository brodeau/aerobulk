/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include "aerobulk.hpp"
#include <iostream>

int main(int argc, char** argv)
{

    std::vector<double> zsst  = {273.15 + 22., 273.15 + 22.};
    std::vector<double> zt_zt = {273.15 + 20., 273.15 + 20.};
    std::vector<double> zq_zt = {0.012, 0.012};
    std::vector<double> zU_zu = {5., 5.};
    std::vector<double> zV_zu = {0., 0.};
    std::vector<double> zslp  = {101000.0, 101000.0};

    std::vector<double> zRsw  = {0., 0.}; // (night)
    std::vector<double> zRlw  = {350., 350.};

    std::vector<double> zQL;
    std::vector<double> zQH;
    std::vector<double> zTau_x;
    std::vector<double> zTau_y;
    std::vector<double> zT_s;

    aerobulk::model(aerobulk::algorithm::COARE, 2, 10, zsst, zt_zt,
            zq_zt, zU_zu, zV_zu, zslp,
            zQL, zQH, zTau_x, zTau_y,
            zRsw, zRlw, zT_s );

    std::cout
        << "\n *********** COARE 3.0 *****************\n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << std::endl
        << " T_s = \t" << zT_s[0] << "\t" << zT_s[1] << std::endl;

    std::vector<double> L = aerobulk::lvap(zsst);

    std::cout << " Evaporation = \t" << -zQL[0]*1e-3/L[0] << "\t" << -zQL[0]*1e-3/L[0] << std::endl;


    aerobulk::model(aerobulk::algorithm::ECMWF, 2, 10, zsst, zt_zt,
            zq_zt, zU_zu, zV_zu, zslp,
            zQL, zQH, zTau_x, zTau_y,
            zRsw, zRlw, zT_s );

    std::cout
        << "\n *********** ECMWF *****************\n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << std::endl
        << " T_s = \t" << zT_s[0] << "\t" << zT_s[1] << std::endl;

    L = aerobulk::lvap(zsst);

    std::cout << " Evaporation = \t" << -zQL[0]*1e-3/L[0] << "\t" << -zQL[0]*1e-3/L[0] << std::endl;

    aerobulk::model(aerobulk::algorithm::NCAR, 2, 10, zsst, zt_zt,
            zq_zt, zU_zu, zV_zu, zslp,
            zQL, zQH, zTau_x, zTau_y);

    std::cout
        << "\n *********** NCAR *****************\n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << std::endl
        << " T_s = \t" << zsst[0] << "\t" << zsst[1] << std::endl;

    L = aerobulk::lvap(zsst);

    std::cout << " Evaporation = \t" << -zQL[0]*1e-3/L[0] << "\t" << -zQL[0]*1e-3/L[0] << std::endl;
}
 

