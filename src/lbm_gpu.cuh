#ifndef LBMPOSTELS_LBM_GPU_CUH
#define LBMPOSTELS_LBM_GPU_CUH

#include <span>

#include "lattice.h"
#include "lbm.h"
#include "types.h"

__global__ void initialize_populations_gpu(std::span<double> f_pop, std::span<double> g_pop, std::span<const double> velocity,
                                std::span<const double> density, std::span<const double> temperature) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i >= N) {
        return;
    }
    auto init_data = get_cell_data(i, velocity, density, temperature);
    for (int q = 0; q < Q; ++q) {
        f_pop[i * Q + q] = f_eq_q(q, init_data);
        g_pop[i * Q + q] = g_eq_q(q, init_data);
    }
}

__global__ void
initialize_moments_gpu(std::span<double> velocity, std::span<double> density, std::span<double> temperature) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i >= N) {
        return;
    }
    for (int dim = 0; dim < 3; ++dim) {
        velocity[i * 3 + dim] = vel_init(i, dim);
    }
    density[i] = density_init(i);
    temperature[i] = temperature_init(i);
}

__global__ void
relax_populations_gpu(std::span<double> f_pop, std::span<double> g_pop, std::span<const double> velocity,
                      std::span<const double> density, std::span<const double> temperature) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i >= Q * N) {
        return;
    }
    auto data = get_cell_data(i / Q, velocity, density, temperature);
    f_pop[i] = (1.0 - relaxation_time_inv) * f_pop[i] + relaxation_time_inv * f_eq_q(i % Q, data);
    g_pop[i] = (1.0 - relaxation_time_g_inv) * g_pop[i] + relaxation_time_g_inv * g_eq_q(i % Q, data);
}

__global__ void stream_population_gpu(std::span<double> f_pop, std::span<double> f_pop_buffer) {
    int cta_idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (cta_idx >= N * Q) {
        return;
    }
    auto [i, j, k, q] = cell_vec_pos(cta_idx);
    f_pop_buffer[cta_idx] = f_pop[cell_pop_idx(i - cqx[q], j - cqy[q], k - cqz[q]) + q];

    if (i == 0 && cqx[q] == 1 || i == Nx - 1 && cqx[q] == -1) {
        f_pop_buffer[cta_idx] = f_pop[cta_idx + (q % 2 ? 1 : -1)];
    }
}

__global__ void update_moments_gpu(std::span<double> velocity, std::span<double> density, std::span<double> temperature,
                                   std::span<const double> f_pop, std::span<const double> g_pop) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx >= N) {
        return;
    }
    update_cell_moment(idx, velocity, density, temperature, f_pop, g_pop);
}

__global__ void apply_force_gpu(std::span<double> velocity) {
    // Gravity?
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i >= N) {
        return;
    }
    for (int dim = 0; dim < 3; ++dim) {
        cell_force(i, dim, velocity[i * 3 + dim]);
    }
}



#endif //LBMPOSTELS_LBM_GPU_CUH
