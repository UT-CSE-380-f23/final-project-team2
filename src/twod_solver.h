#pragma once
#include "fd_solver.h"

constexpr double scaling_constant_2d = 1.0 / (8.0 * M_PI * M_PI); ;
class TwoDSolver:public FDSolver{
    public:
        TwoDSolver();
        TwoDSolver(const size_t& num_nodes, const bool& solver_method);
        TwoDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order);
        TwoDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order,const bool& verify, const bool& debug);
        void construct_matrix();
};