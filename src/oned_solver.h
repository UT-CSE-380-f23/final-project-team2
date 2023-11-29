#pragma once
#include "fd_solver.h"

constexpr double scaling_constant = 1.0 / (4 * M_PI * M_PI);
class OneDSolver:public FDSolver{
    public:
        OneDSolver();
        OneDSolver(const size_t& num_nodes, const bool& solver_method);
        OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order);
        OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order,const bool& verify, const bool& debug);
        void construct_matrix();
};
