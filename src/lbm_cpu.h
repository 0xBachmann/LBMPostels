#ifndef LBMPOSTELS_LBM_CPU_H
#define LBMPOSTELS_LBM_CPU_H

#include "constants.h"
#include "types.h"
#include "lbm.h"
#include "vtk.h"

void initialize_populations(PopulationT& f_pop, PopulationT& g_pop, const VelocityT& velocity, const ScalarT& density, const ScalarT& temperature) {
    for(int i = 0; i < N; ++i) {
        auto init_data = get_cell_data(i, velocity, density, temperature);
        for(int q = 0; q < Q; ++q) {
            f_pop[i*Q + q] = f_eq_q(q, init_data);
            g_pop[i*Q + q] = g_eq_q(q, init_data);
        }
    }
}

void initialize_moments(VelocityT& velocity, ScalarT& density, ScalarT& temperature) {
    for(int i = 0; i < N; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            velocity[i * 3 + dim] = vel_init(i, dim);
        }
        density[i] = density_init(i);
        temperature[i] = temperature_init(i);
    }
}

void relax_populations(PopulationT& f_pop, PopulationT& g_pop, const VelocityT& velocity, const ScalarT& density, const ScalarT& temperature) {
    for(int i = 0; i < N; ++i) {
        auto data = get_cell_data(i, velocity, density, temperature);
        for(int q = 0; q < Q; ++q) {
            f_pop[i*Q + q] = (1.0 - relaxation_time_inv) * f_pop[i*Q + q] + relaxation_time_inv * f_eq_q(q, data);
            g_pop[i*Q + q] = (1.0 - relaxation_time_g_inv) * g_pop[i*Q + q] + relaxation_time_g_inv * g_eq_q(q, data);
        }
    }
}

void stream_population(PopulationT& f_pop, PopulationT& f_pop_buffer) {
    // Inner cells
    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                const int current_cell_idx = cell_pop_idx(i, j, k);
                for(int q = 0; q < Q; ++q) {
                    f_pop_buffer[current_cell_idx + q] = f_pop[cell_pop_idx(i - cqx[q], j - cqy[q], k - cqz[q]) + q];
                }
            }
        }
    }

    // Boundaries
    for(int j = 1; j < Ny - 1 ; ++j){
        for(int k = 0; k < Nz; ++k) {
            // Left x-
            int i = 0;
            apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqx[q] == 1;});

            // Right x+
            i = Nx - 1;
            apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqx[q] == -1;});
        }
    }
    for(int i = 1; i < Nx - 1; ++i) {
        for(int k = 0; k < Nz; ++k) {
            // Bottom: y-
            int j = 0;
            // apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqy[q] == 1;});

            // Top: y+
            j = Ny - 1;
            // apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqy[q] == -1;});
        }

        // for(int j = 1; j < Ny; ++j) {
        //     // Back: z-
        //     int k = 0;
        //     apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqz[q] == 1;});

        //     // Front: z+
        //     k = Nz - 1;
        //     apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqz[q] == -1;});
        // }
    }

    std::swap(f_pop, f_pop_buffer);


    // // Corners
    // // Left Bottom Back
    // apply_BB_boundary(f_pop, f_pop_buffer, 0, 0, 0, [](int q){return cqx[q] == 1 || ;});


}

void update_moments(VelocityT& velocity, ScalarT& density, ScalarT& temperature, const PopulationT& f_pop, const PopulationT& g_pop) {
    for(int i = 0; i < N; ++i) {
        update_cell_moment(i, velocity, density, temperature, f_pop, g_pop);
    }
}

void apply_force(VelocityT& velocity) {
    // Gravity?
    for(int i = 0; i < N; ++i) {
        for (int dim = 0; dim < 3; ++dim) {
            cell_force(i, dim, velocity[i*3 + dim]);
        }
    }
}

#endif //LBMPOSTELS_LBM_CPU_H
