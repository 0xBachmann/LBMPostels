#ifndef LBMPOSTELS_CONSTANTS_H
#define LBMPOSTELS_CONSTANTS_H

#include <cstddef>
#include <array>

#include "cuda_utils.h"

constexpr CUDA_CONST int output_frequency = 10;

constexpr CUDA_CONST std::size_t Nx = 200;
constexpr CUDA_CONST std::size_t Ny = 200;
constexpr CUDA_CONST std::size_t Nz = 1;
constexpr CUDA_CONST std::size_t N = Nx*Ny*Nz;
constexpr CUDA_CONST std::size_t Q = 19;
constexpr CUDA_CONST std::size_t Npop = N*Q;


// TODO: switch cix dependent on Q
constexpr CUDA_CONST std::array<double, Q> cqx{0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
constexpr CUDA_CONST std::array<double, Q> cqy{0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
constexpr CUDA_CONST std::array<double, Q> cqz{0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
constexpr CUDA_CONST std::array<double, Q> wq{1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};

constexpr CUDA_CONST double dt = 1.;
constexpr CUDA_CONST double dx = dt;


// Speed of sound (in lattice units)
constexpr CUDA_CONST double cs_2 = 3.; // cs^2 = 1./3. * (dx**2/dt**2)
constexpr CUDA_CONST double cs_4 = 9.;

constexpr CUDA_CONST double viscosity = 5e-2;
constexpr CUDA_CONST double D = 1e-2;
constexpr CUDA_CONST double U = 0.2;
constexpr CUDA_CONST double R = Nx/5;

constexpr CUDA_CONST double Re = R*U/viscosity;

constexpr CUDA_CONST double relaxation_time = viscosity * cs_2 + 0.5;
constexpr CUDA_CONST double relaxation_time_g = D * cs_2 + 0.5;
constexpr CUDA_CONST double relaxation_time_inv = 1.0 / relaxation_time;
constexpr CUDA_CONST double relaxation_time_g_inv = 1.0 / relaxation_time_g;

constexpr CUDA_CONST bool incompressible = false;
constexpr CUDA_CONST double rho0 = 1.0;

#endif //LBMPOSTELS_CONSTANTS_H
