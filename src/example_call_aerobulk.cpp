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

    aerobulk::model(std::string("coare"), 2, 10, zsst, zt_zt,
            zq_zt, zU_zu, zV_zu, zslp,
            zQL, zQH, zTau_x, zTau_y,
            10, zRsw, zRlw );

    std::cout
        << " QH=" << zQH[0] << " " << zQH[1] << std::endl
        << " QL=" << zQL[0] << " " << zQL[1] << std::endl
        << " Tau_x=" << zTau_x[0] << " " << zTau_x[1] << std::endl
        << " Tau_y=" << zTau_y[0] << " " << zTau_y[1] << std::endl;
}
 

