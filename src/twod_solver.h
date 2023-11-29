#pragma once
#include "fd_solver.h"
/*
** TwoDSolver class is a child class of FDSolver containing only the matrix construction in 2D 2nd and 4th order
** and constructors that help initialize the FD Solver class
*/

inline constexpr double scaling_constant_2d = 1.0 / (8.0 * M_PI * M_PI); ;
class TwoDSolver:public FDSolver{
    public:
        TwoDSolver();
        TwoDSolver(const size_t& num_nodes, const bool& solver_method);
        TwoDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order);
        TwoDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order,const bool& verify, const bool& debug);
        void construct_matrix();
};