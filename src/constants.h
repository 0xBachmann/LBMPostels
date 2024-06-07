#ifndef LBMPOSTELS_CONSTANTS_H
#define LBMPOSTELS_CONSTANTS_H

#include <cstddef>
#include <array>

constexpr int output_frequency = 10;

constexpr std::size_t Nx = 40;
constexpr std::size_t Ny = 70;
constexpr std::size_t Nz = 1;
constexpr std::size_t N = Nx*Ny*Nz;
constexpr std::size_t Q = 19;
constexpr std::size_t Npop = N*Q;


// TODO: switch cix dependent on Q
constexpr std::array<double, Q> cqx{0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
constexpr std::array<double, Q> cqy{0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
constexpr std::array<double, Q> cqz{0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
constexpr std::array<double, Q> wq{1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};

constexpr double dt = 1.;
constexpr double dx = dt;


// Speed of sound (in lattice units)
constexpr double cs_2 = 3.; // cs^2 = 1./3. * (dx**2/dt**2)
constexpr double cs_4 = 9.;

constexpr double viscosity = 5e-2;
constexpr double D = 1e-2;
constexpr double U = 0.2;
constexpr double R = Nx/5;

constexpr double Re = R*U/viscosity;

constexpr double relaxation_time = viscosity * cs_2 + 0.5;
constexpr double relaxation_time_g = D * cs_2 + 0.5;
constexpr double relaxation_time_inv = 1.0 / relaxation_time;
constexpr double relaxation_time_g_inv = 1.0 / relaxation_time_g;

constexpr bool incompressible = false;
constexpr double rho0 = 1.0;

#endif //LBMPOSTELS_CONSTANTS_H
