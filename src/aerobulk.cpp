#include "aerobulk.hpp"

extern "C"
{
	void aerobulk_cxx_all( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const int *, const double *, const double *,
			const int *, const int *);

	void aerobulk_cxx_rad( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const double *, const double *,
			const int *, const int *);

	void aerobulk_cxx_niter( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const int *,
			const int *, const int *);

	void aerobulk_cxx_noopt( const char *, const double *, const double *, const double *, const double *,
			const double *, const double *, const double *, const double *,
			double *, double *, double *, double *,
			const int *, const int *);
}

void aerobulk::check_inputs(int count, ...)
{
    va_list ap;

    va_start(ap, count); // Requires the last fixed parameter (to get the address)
    int m = va_arg(ap, int); // Get the first size - the one we compare everyone else's with

    for (int i=1; i<count; i++) // Start from 1 because we already called var_arg once
        assert( m == va_arg(ap, int) ); // Increments ap to the next argument

    va_end(ap);
}

void aerobulk::model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
		const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
		std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
		int Niter, const std::vector<double> &rad_sw, const std::vector<double> &rad_lw)
{
	int l = calgo.size();
	int m = sst.size();

	aerobulk::check_inputs(12, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size(), rad_sw.size(), rad_lw.size(),
            QL.size(), QH.size(), Tau_x.size(), Tau_y.size());

	aerobulk_cxx_all(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
            &QL[0], &QH[0], &Tau_x[0], &Tau_y[0],
			&Niter, &rad_sw[0], &rad_lw[0],
			&l, &m);
}

void aerobulk::model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
		const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
		std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
		const std::vector<double> &rad_sw, const std::vector<double> &rad_lw)
{
	int l = calgo.size();
	int m = sst.size();

	aerobulk::check_inputs(8, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size(), rad_sw.size(), rad_lw.size());

	aerobulk_cxx_rad(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
            &QL[0], &QH[0], &Tau_x[0], &Tau_y[0],
			&rad_sw[0], &rad_lw[0],
			&l, &m);
}

void aerobulk::model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
		const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
		std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
		const int Niter)
{
	int l = calgo.size();
	int m = sst.size();

	aerobulk::check_inputs(6, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size());

    aerobulk_cxx_niter(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
            &QL[0], &QH[0], &Tau_x[0], &Tau_y[0],
            &Niter,
            &l, &m);
}

void aerobulk::model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
		const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
		std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y)
{
	int l = calgo.size();
	int m = sst.size();

	aerobulk::check_inputs(6, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), V_zu.size(), slp.size());

	aerobulk_cxx_noopt(calgo.c_str(), &zt, &zu, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &V_zu[0], &slp[0],
            &QL[0], &QH[0], &Tau_x[0], &Tau_y[0],
			&l, &m);
}

