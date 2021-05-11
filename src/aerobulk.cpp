/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include "aerobulk.hpp"

extern "C"
{
	void aerobulk_cxx_skin( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const double *, const double *, double *,
			const int *, const int *, const int *);

	void aerobulk_cxx_no_skin( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const int *, const int *, const int *);

    void l_vap_cxx( const double *, int *, double *);
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
        case algorithm::COARE3p6:
            return_value = std::string("coare3p6");
            break;
        case algorithm::NCAR:
            return_value = std::string("ncar");
            break;
        case algorithm::ECMWF:
            return_value = std::string("ecmwf");
            break;
        default:
            return_value = "unknown";
    }
    return return_value;
}

// Check the size of input data
int aerobulk::check_sizes(int count, ...)
{
    va_list ap;

    va_start(ap, count); // Requires the last fixed parameter (to get the address)
    int size = va_arg(ap, int); // Get the first size - the one we compare everyone else's with

    for (int i=1; i<count; i++) // Start from 1 because we already called var_arg once
        assert( size == va_arg(ap, int) ); // Increments ap to the next argument

    va_end(ap);

    return size;
}

// Interface for l_vap
std::vector<double> aerobulk::l_vap(const std::vector<double> &sst)
{
    // Prepp
    int m = sst.size();
    std::vector<double> l_vap_out(m);

    // The actual function call
    l_vap_cxx(&sst[0], &m, &l_vap_out[0]);

    return l_vap_out;
}

// Interface to aerobulk_model with rad_sw and rad_lw as inputs and T_s as output
void aerobulk::model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
    const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
    std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
    const std::vector<double> &rad_sw, const std::vector<double> &rad_lw, std::vector<double> &T_s,
    const int Niter)
{
    // Algorithm string and size
    std::string calgo = aerobulk::algorithm_to_string(algo);
	int l = calgo.size();

    // Check the length of the inputs and record it
	int m = aerobulk::check_sizes(8, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size(), rad_sw.size(), rad_lw.size());

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Tau_x.resize(m);
    Tau_y.resize(m);
    T_s.resize(m);

    // The actual function call - we need to sent the adresses/pointer because it's a C interface to a Fortran routine
    aerobulk_cxx_skin(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
        &QL[0], &QH[0], &Tau_x[0], &Tau_y[0],
        &rad_sw[0], &rad_lw[0], &T_s[0], &Niter,
        &l, &m);
}

// Interface to aerobulk_model without rad_sw, rad_lw, and T_s
void aerobulk::model(algorithm algo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
    const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
    std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
    const int Niter)
{
    // Algorithm string and size
    std::string calgo = aerobulk::algorithm_to_string(algo);
	int l = calgo.size();

    // Check the length of the inputs and record it
	int m = aerobulk::check_sizes(6, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size());

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Tau_x.resize(m);
    Tau_y.resize(m);

    // The actual function call - we need to sent the adresses/pointer because it's a C interface to a Fortran routine
    aerobulk_cxx_no_skin(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
            &QL[0], &QH[0], &Tau_x[0], &Tau_y[0], &Niter,
            &l, &m);
}

