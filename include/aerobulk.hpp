#include <vector>
#include <string>
#include <cassert>
#include <cstdarg>


namespace aerobulk
{

    // To check the size of the inputs
    int check_sizes(int count, ...);

    // Interface for all options activated
    void model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            int Niter, const std::vector<double> &rad_sw, const std::vector<double> &rad_lw);

    // Interface for only using optional incoming radiation
    void model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            const std::vector<double> &rad_sw, const std::vector<double> &rad_lw);

    // Interface for only using optional number of iterations
    void model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y,
            int Niter);

    // Interface for only using none of the optional inputs
    void model(std::string calgo, double zt, double zu, const std::vector<double> &sst, const std::vector<double> &t_zt,
            const std::vector<double> &q_zt, const std::vector<double> &U_zu, const std::vector<double> &V_zu, const std::vector<double> &slp,
            std::vector<double> &QL, std::vector<double> &QH, std::vector<double> &Tau_x, std::vector<double> &Tau_y);

}
