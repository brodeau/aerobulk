/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include "aerobulk.hpp"

extern "C"
{
    void aerobulk_cxx_skin(    const int *, const int *, const char *, const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, const double *,
                               double *, double *, double *, double *, double *,
                               const int *, const bool *, const double *, const double *, double *,
                               const int *, const int * );

    // Same as aerobulk_cxx_skin, plus the UTC time (seconds since 00:00 UTC) and the longitude,
    // needed by the warm-layer scheme of COARE
    void aerobulk_cxx_skin_wl( const int *, const int *, const char *, const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, const double *,
                               double *, double *, double *, double *, double *,
                               const int *, const bool *, const double *, const double *, double *,
                               const int *, const double *,
                               const int *, const int * );

    void aerobulk_cxx_no_skin( const int *, const int *, const char *, const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, const double *,
                               double *, double *, double *, double *, double *,
                               const int *,
                               const int *, const int * );

}

// Convert the C++ algorithm enumeration to a string
std::string aerobulk::algorithm_to_string(aerobulk::algorithm algo)
{
    std::string return_value;
    switch (algo)
    {
    case algorithm::OTHER:
        return_value = std::string("other");
        break;
    case algorithm::COARE3p0:
        return_value = std::string("coare3p0");
        break;
    case algorithm::COARE3p6:
        return_value = std::string("coare3p6");
        break;
    case algorithm::NCAR:
        return_value = std::string("ncar");
        break;
    case algorithm::ECMWF:
        return_value = std::string("ecmwf");
        break;
    case algorithm::ANDREAS:
        return_value = std::string("andreas");
        break;
    default:
        return_value = "unknown";
    }
    return return_value;
}

// Check that all input vectors have the same size, and return that size.
// We use a real check that throws an exception, not `assert`, because `assert` is removed
// when NDEBUG is defined (release builds), and a size mismatch would then make Fortran read out of bounds.
// The sizes are taken as `size_t` (the type of std::vector::size()), so there is no type mismatch
// (the old variadic version read them as `int` with va_arg, which is undefined behaviour).
int aerobulk::check_sizes(std::initializer_list<std::size_t> sizes)
{
    const std::size_t size = *sizes.begin(); // the size we compare all the others with

    for (const std::size_t s : sizes)
        if (s != size)
            throw std::invalid_argument("aerobulk::check_sizes => input vectors do not all have the same size!");

    return static_cast<int>(size);
}

// Interface for l_vap
/*
  std::vector<double> aerobulk::l_vap(const std::vector<double> &sst)
  {
  // Prepp
  int m = sst.size();
  std::vector<double> l_vap_out(m);

  // The actual function call
  l_vap_cxx(&sst[0], &m, &l_vap_out[0]);

  return l_vap_out;
  }
*/

// Interface to aerobulk_model with rad_sw and rad_lw as inputs and T_s as output
void aerobulk::model(const int jt, const int Nt, algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
                     const std::vector<double> &hum_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
                     std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y, std::vector<double> &Evap,
                     const int Niter, const bool l_use_skin, const std::vector<double> &rad_sw, const std::vector<double> &rad_lw, std::vector<double> &T_s)

{
    // Algorithm string and size
    std::string calgo = aerobulk::algorithm_to_string(algo);
    int l = calgo.size();

    // Check the length of the inputs and record it
    int m = aerobulk::check_sizes({sst.size(), t_zt.size(), hum_zt.size(), U_zu.size(), V_zu.size(), slp.size(), rad_sw.size(), rad_lw.size()});

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Tau_x.resize(m);
    Tau_y.resize(m);
    Evap.resize(m);
    T_s.resize(m);

    // The actual function call - we need to sent the adresses/pointer because it's a C interface to a Fortran routine
    aerobulk_cxx_skin( &jt, &Nt, calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &hum_zt[0], &U_zu[0], &V_zu[0], &slp[0],
                       &QL[0], &QH[0], &Tau_x[0], &Tau_y[0], &Evap[0],
                       &Niter, &l_use_skin, &rad_sw[0], &rad_lw[0], &T_s[0],
                       &l, &m );
}


// Interface to aerobulk_model with rad_sw and rad_lw as inputs and T_s as output,
// plus the UTC time and the longitude, needed by the warm-layer scheme of COARE
void aerobulk::model(const int jt, const int Nt, algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
                     const std::vector<double> &hum_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
                     std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y, std::vector<double> &Evap,
                     const int Niter, const bool l_use_skin, const std::vector<double> &rad_sw, const std::vector<double> &rad_lw, std::vector<double> &T_s,
                     const int isecday_utc, const std::vector<double> &plong)
{
    // Algorithm string and size
    std::string calgo = aerobulk::algorithm_to_string(algo);
    int l = calgo.size();

    // Check the length of the inputs (including the longitude array) and record it
    int m = aerobulk::check_sizes({sst.size(), t_zt.size(), hum_zt.size(), U_zu.size(), V_zu.size(), slp.size(), rad_sw.size(), rad_lw.size(), plong.size()});

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Tau_x.resize(m);
    Tau_y.resize(m);
    Evap.resize(m);
    T_s.resize(m);

    // The actual function call - we need to send the adresses/pointers because it's a C interface to a Fortran routine
    aerobulk_cxx_skin_wl( &jt, &Nt, calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &hum_zt[0], &U_zu[0], &V_zu[0], &slp[0],
                          &QL[0], &QH[0], &Tau_x[0], &Tau_y[0], &Evap[0],
                          &Niter, &l_use_skin, &rad_sw[0], &rad_lw[0], &T_s[0],
                          &isecday_utc, &plong[0],
                          &l, &m );
}




// Interface to aerobulk_model without rad_sw, rad_lw, and T_s
void aerobulk::model( const int jt, const int Nt, algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
                      const std::vector<double> &hum_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp ,
                      std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y, std::vector<double> &Evap,
                      const int Niter)
{
    // Algorithm string and size
    std::string calgo = aerobulk::algorithm_to_string(algo);
    int l = calgo.size();

    // Check the length of the inputs and record it
    int m = aerobulk::check_sizes({sst.size(), t_zt.size(), hum_zt.size(), U_zu.size(), V_zu.size(), slp.size()});

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Tau_x.resize(m);
    Tau_y.resize(m);
    Evap.resize(m);

    // The actual function call - we need to sent the adresses/pointer because it's a C interface to a Fortran routine
    aerobulk_cxx_no_skin( &jt, &Nt, calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &hum_zt[0], &U_zu[0], &V_zu[0], &slp[0],
                          &QL[0], &QH[0], &Tau_x[0], &Tau_y[0], &Evap[0],
                          &Niter, &l, &m );
}
