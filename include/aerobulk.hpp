/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#ifndef __aerobulk_HPP
#define __aerobulk_HPP 1

#include <vector>
#include <string>
#include <cassert>
#include <cstdarg>

namespace aerobulk
{
    enum class algorithm
    {
        COARE   = 0,
        COARE35 = 1,
        NCAR    = 2,
        ECMWF   = 3
    };

    std::string algorithm_to_string(algorithm algo);

    // To check the size of the inputs
    int check_sizes(int count, ...);

    // Interface for all options activated
    void model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            int Niter, const std::vector<double> &rad_sw, const std::vector<double> &rad_lw);

    // Interface for only using optional incoming radiation
    void model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            const std::vector<double> &rad_sw, const std::vector<double> &rad_lw);

    // Interface for only using optional number of iterations
    void model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            int Niter);

    // Interface for only using none of the optional inputs
    void model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y);

}

#endif
