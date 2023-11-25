#ifndef FD_SOLVER_H
#define FD_SOLVER_H
#include "fd_solver.h"

inline constexpr double scaling_constant = 1.0 / (4 * M_PI * M_PI);
class OneDSolver:public FDSolver{
    public:
        OneDSolver();
        OneDSolver(const size_t& num_nodes, const bool& solver_method);
        OneDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order);
        
        void construct_matrix();
};
#endif