#ifndef LBMPOSTELS_LBM_H
#define LBMPOSTELS_LBM_H

#include "types.h"

#include "types.h"
#include "constants.h"
#include "lattice.h"

CUDA_HOST_DEVICE double
f_eq_q_compressible(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    const double uc = ux * cqx[q] + uy * cqy[q] + uz * cqz[q];
    const double uc2 = uc * uc;
    return wq[q] * density * (1. + uc * cs_2 + 0.5 * uc2 * cs_4 - 0.5 * u2 * cs_2);
}

CUDA_HOST_DEVICE double
g_eq_q(int q, const double ux, const double uy, const double uz, const double u2, const double T) {
    return f_eq_q_compressible(q, ux, uy, uz, u2, T);
}

CUDA_HOST_DEVICE double
f_eq_q_incompressible(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    const double uc = ux * cqx[q] + uy * cqy[q] + uz * cqz[q];
    const double uc2 = uc * uc;
    return wq[q] * density + wq[q] * rho0 * (uc * cs_2 + 0.5 * uc2 * cs_4 - 0.5 * u2 * cs_2);
}

CUDA_HOST_DEVICE double
f_eq_q(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    if constexpr (incompressible) {
        return f_eq_q_incompressible(q, ux, uy, uz, u2, density);
    } else {
        return f_eq_q_compressible(q, ux, uy, uz, u2, density);
    }
}

struct cell_data {
    double ux;
    double uy;
    double uz;
    double u2;
    double cell_density;
    double cell_temperature;
};

CUDA_HOST_DEVICE double f_eq_q(int q, const cell_data &init_data) {
    return f_eq_q(q, init_data.ux, init_data.uy, init_data.uz, init_data.u2, init_data.cell_density);
}

CUDA_HOST_DEVICE double g_eq_q(int q, const cell_data &init_data) {
    return g_eq_q(q, init_data.ux, init_data.uy, init_data.uz, init_data.u2, init_data.cell_temperature);
}

CUDA_HOST_DEVICE cell_data
get_cell_data(int cell_idx, std::span<const double> velocity, std::span<const double> density,
              std::span<const double> temperature) {
    return {
            .ux = velocity[cell_idx * 3],
            .uy = velocity[cell_idx * 3 + 1],
            .uz = velocity[cell_idx * 3 + 2],
            .u2 = velocity[cell_idx * 3] * velocity[cell_idx * 3] +
                  velocity[cell_idx * 3 + 1] * velocity[cell_idx * 3 + 1] +
                  velocity[cell_idx * 3 + 2] * velocity[cell_idx * 3 + 2],
            .cell_density = density[cell_idx],
            .cell_temperature = temperature[cell_idx]
    };
}

CUDA_HOST_DEVICE double vel_init(int idx, int dim) {
    auto [i, j, k] = cell_pos(idx);

    if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 3) * (j - Ny / 3) < (R) * (R)) {
        return 0;
    } else {
        return dim == 1 ? U : 0;
    }
}

CUDA_HOST_DEVICE double density_init(int idx) {
    return 1.0;
}

CUDA_HOST_DEVICE double temperature_init(int idx) {
    auto [i, j, k] = cell_pos(idx);

    if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 3) * (j - Ny / 3) < (R) * (R)) {
        return 1.0;
    } else {
        return 0.0;
    }
}

CUDA_HOST_DEVICE void
update_cell_moment(int cell_idx, std::span<double> velocity, std::span<double> density, std::span<double> temperature, std::span<const double> f_pop,
                   std::span<const double> g_pop) {
    double cell_density = 0.;
    double ux = 0.;
    double uy = 0.;
    double uz = 0.;
    double cell_temperature = 0.;

    for (int q = 0; q < Q; ++q) {
        cell_density += f_pop[cell_idx * Q + q];
        cell_temperature += g_pop[cell_idx * Q + q];
        ux += cqx[q] * f_pop[cell_idx * Q + q];
        uy += cqy[q] * f_pop[cell_idx * Q + q];
        uz += cqz[q] * f_pop[cell_idx * Q + q];
    }

    double cell_density_inv = 1. / cell_density;
    ux *= cell_density_inv;
    uy *= cell_density_inv;
    uz *= cell_density_inv;

    velocity[cell_idx * 3] = ux;
    velocity[cell_idx * 3 + 1] = uy;
    velocity[cell_idx * 3 + 2] = uz;

    density[cell_idx] = cell_density;
    temperature[cell_idx] = cell_temperature;
}

CUDA_HOST_DEVICE void cell_force(int idx, int dim, double &vel) {
    auto [i, j, k] = cell_pos(idx);
    if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 3) * (j - Ny / 3) < (R) * (R)) {
        vel = 0;
    }
}

// BB Boundary
// If q is odd we have to add one to it to find the flipped vector index
// If even we remove one
template<typename P>
CUDA_HOST_DEVICE void apply_BB_boundary(std::span<double> f_pop, std::span<double> f_pop_buffer, int i, int j, int k, P&& predicate) {
    const double current_cell_idx = cell_pop_idx(i, j, k);
    for(int q = 0; q < Q; ++q) {
        if(predicate(q)) {
            f_pop_buffer[current_cell_idx + q] = f_pop[current_cell_idx + (q%2 ? q + 1 : q - 1)];
        } else {
            f_pop_buffer[current_cell_idx + q] = f_pop[cell_pop_idx(i - cqx[q], j - cqy[q], k - cqz[q]) + q];
        }
    }
}

#endif //LBMPOSTELS_LBM_H
