#pragma once
#include "fd_solver.h"

/*
** OneDSolver class is a child class of FDSolver containing only the matrix construction in 1D 2nd and 4th order
** and constructors that help initialize the FD Solver class
*/
inline constexpr double scaling_constant = 1.0 / (4 * M_PI * M_PI);
class OneDSolver:public FDSolver{
    public:
        OneDSolver();
        OneDSolver(const size_t& num_nodes, const bool& solver_method);
        OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order);
        OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order,const bool& verify, const bool& debug);
        void construct_matrix();
};
