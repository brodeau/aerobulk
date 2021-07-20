/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include "aerobulk.hpp"
#include <iostream>

int main(int argc, char** argv)
{

    double zt =  2.; // height of measurement for air temperature and humidity [m]
    double zu = 10.; // height of measurement for wind speed [m]

    int nbiter = 8 ; // number of iterations...

    double to_mm_p_day =  3600.*24.;

    std::vector<double> zsst  = {273.15 + 22., 273.15 + 22.}; // sea surface temperature [K]
    std::vector<double> zt_zt = {273.15 + 20., 273.15 + 25.}; // air absolute temperature at zt [K] ## second case is stable ABL as t_air > SST (25>22)!
    std::vector<double> zq_zt = {0.012, 0.012};               // air specific humidity at zt [g/kg]
    std::vector<double> zU_zu = { 4.,  4.};                   // u-wind speed component at zu [m/s]
    std::vector<double> zV_zu = { 9.,  9.};                   // v-wind speed component at zu [m/s]
    std::vector<double> zslp  = {101000.0, 101000.0};         // sea-level atmospheric pressure [Pa]

    std::vector<double> zRsw  = {0., 0.};       // downwelling shortwave (solar)     radiation [W/m^2] ## night!
    std::vector<double> zRlw  = {350., 350.};   // downwelling longwave  (infra-red) radiation [W/m^2]

    // AeroBulk output:
    std::vector<double> zQL;
    std::vector<double> zQH;
    std::vector<double> zTau_x;
    std::vector<double> zTau_y;
    std::vector<double> zE;
    std::vector<double> zT_s;


    /* ############################################################################################## */

    std::cout << "\n\n *********** COARE 3.0 parameterization *****************" << std::endl;

    aerobulk::model( 1, 1, aerobulk::algorithm::COARE3p0, zt, zu, zsst, zt_zt,
                     zq_zt, zU_zu, zV_zu, zslp,
                     zQL, zQH, zTau_x, zTau_y, zE, 
                     nbiter, true, zRsw, zRlw, zT_s );

    std::cout
        << "  ==> results for COARE 3.0: \n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << "    [W/m^2]" << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << "    [W/m^2]" << std::endl
        << " Evap = \t" << to_mm_p_day*zE[0] << "\t" << to_mm_p_day*zE[1] << "    [mm/day]" << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << "    [N/m^2]" << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << "    [N/m^2]" << std::endl
        << " T_s = \t" << zT_s[0] << "\t" << zT_s[1] << "    [K]" << std::endl;

    // std::vector<double> L = aerobulk::l_vap(zsst);
    //std::cout << " Evaporation = \t" << -zQL[0]*1e-3/L[0] << "\t" << -zQL[1]*1e-3/L[1] << std::endl;



    /* ############################################################################################## */

    std::cout << "\n\n *********** COARE 3.6 parameterization *****************" << std::endl;

    aerobulk::model( 1, 1, aerobulk::algorithm::COARE3p6, zt, zu, zsst, zt_zt,
                     zq_zt, zU_zu, zV_zu, zslp,
                     zQL, zQH, zTau_x, zTau_y, zE,
                     nbiter, true, zRsw, zRlw, zT_s );

    std::cout
        << "  ==> results for COARE 3.6: \n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << "    [W/m^2]" << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << "    [W/m^2]" << std::endl
        << " Evap = \t" << to_mm_p_day*zE[0] << "\t" << to_mm_p_day*zE[1] << "    [mm/day]" << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << "    [N/m^2]" << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << "    [N/m^2]" << std::endl
        << " T_s = \t" << zT_s[0] << "\t" << zT_s[1] << "    [K]" << std::endl;



    /* ############################################################################################## */

    std::cout << "\n\n *********** ECMWF parameterization *****************\n" << std::endl;

    aerobulk::model( 1, 1, aerobulk::algorithm::ECMWF, zt, zu, zsst, zt_zt,
                     zq_zt, zU_zu, zV_zu, zslp,
                     zQL, zQH, zTau_x, zTau_y, zE,
                     nbiter, true, zRsw, zRlw, zT_s );

    std::cout
        << "\n  ==> results for ECMWF: \n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << "    [W/m^2]" << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << "    [W/m^2]" << std::endl
        << " Evap = \t" << to_mm_p_day*zE[0] << "\t" << to_mm_p_day*zE[1] << "    [mm/day]" << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << "    [N/m^2]" << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << "    [N/m^2]" << std::endl
        << " T_s = \t" << zT_s[0] << "\t" << zT_s[1] << "    [K]" << std::endl;


    /* ############################################################################################## */

    std::cout << "\n\n *********** NCAR parameterization *****************\n" << std::endl;

    aerobulk::model( 1, 1, aerobulk::algorithm::NCAR, zt, zu, zsst, zt_zt,
                     zq_zt, zU_zu, zV_zu, zslp,
                     zQL, zQH, zTau_x, zTau_y, zE,
                     nbiter );

    std::cout
        << "\n  ==> results for NCAR: \n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << "    [W/m^2]" << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << "    [W/m^2]" << std::endl
        << " Evap = \t" << to_mm_p_day*zE[0] << "\t" << to_mm_p_day*zE[1] << "    [mm/day]" << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << "    [N/m^2]" << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << "    [N/m^2]" << std::endl
        << " T_s = \t" << zsst[0] << "\t" << zsst[1] << "    [K]" << std::endl;


    /* ############################################################################################## */

    std::cout << "\n\n *********** ANDREAS parameterization *****************\n" << std::endl;

    aerobulk::model( 1, 1, aerobulk::algorithm::ANDREAS, zt, zu, zsst, zt_zt,
                     zq_zt, zU_zu, zV_zu, zslp,
                     zQL, zQH, zTau_x, zTau_y, zE,
                     nbiter);

    std::cout
        << "\n  ==> results for ANDREAS: \n"
        << " QH = \t" << zQH[0] << "\t" << zQH[1] << "    [W/m^2]" << std::endl
        << " QL = \t" << zQL[0] << "\t" << zQL[1] << "    [W/m^2]" << std::endl
        << " Evap = \t" << to_mm_p_day*zE[0] << "\t" << to_mm_p_day*zE[1] << "    [mm/day]" << std::endl
        << " Tau_x = \t" << zTau_x[0] << "\t" << zTau_x[1] << "    [N/m^2]" << std::endl
        << " Tau_y = \t" << zTau_y[0] << "\t" << zTau_y[1] << "    [N/m^2]" << std::endl
        << " T_s = \t" << zsst[0] << "\t" << zsst[1] << "    [K]" << std::endl;



}
