#ifndef LBMPOSTELS_LBM_H
#define LBMPOSTELS_LBM_H

#include "types.h"

#include "types.h"
#include "constants.h"
#include "lattice.h"

double
f_eq_q_compressible(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    const double uc = ux * cqx[q] + uy * cqy[q] + uz * cqz[q];
    const double uc2 = uc * uc;
    return wq[q] * density * (1. + uc * cs_2 + 0.5 * uc2 * cs_4 - 0.5 * u2 * cs_2);
}

double g_eq_q(int q, const double ux, const double uy, const double uz, const double u2, const double T) {
    return f_eq_q_compressible(q, ux, uy, uz, u2, T);
}

double
f_eq_q_incompressible(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    const double uc = ux * cqx[q] + uy * cqy[q] + uz * cqz[q];
    const double uc2 = uc * uc;
    return wq[q] * density + wq[q] * rho0 * (uc * cs_2 + 0.5 * uc2 * cs_4 - 0.5 * u2 * cs_2);
}

double f_eq_q(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    if constexpr (incompressible) {
        return f_eq_q_incompressible(q, ux, uy, uz, u2, density);
    } else {
        return f_eq_q_compressible(q, ux, uy, uz, u2, density);
    }
}

struct cell_init_data {
    double ux;
    double uy;
    double uz;
    double u2;
    double cell_density;
    double cell_temperature;
};

double f_eq_q(int q, const cell_init_data &init_data) {
    return f_eq_q(q, init_data.ux, init_data.uy, init_data.uz, init_data.u2, init_data.cell_density);
}

double g_eq_q(int q, const cell_init_data &init_data) {
    return g_eq_q(q, init_data.ux, init_data.uy, init_data.uz, init_data.u2, init_data.cell_temperature);
}

cell_init_data
get_cell_init_data(int cell_idx, const VelocityT &velocity, const ScalarT &density, const ScalarT &temperature) {
    return {
            .ux = velocity[cell_idx * 3],
            .uy = velocity[cell_idx * 3 + 1],
            .uz = velocity[cell_idx * 3 + 2],
            .u2 = velocity[cell_idx] * velocity[cell_idx * 3] +
                  velocity[cell_idx * 3 + 1] * velocity[cell_idx * 3 + 1] +
                  velocity[cell_idx * 3 + 2] * velocity[cell_idx * 3 + 2],
            .cell_density = density[cell_idx],
            .cell_temperature = temperature[cell_idx]
    };
}

double vel_init(int idx, int dim) {
    auto [i, j, k] = cell_pos(idx);

    if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 3) * (j - Ny / 3) < (R) * (R)) {
        return 0;
    } else {
        return dim == 1 ? U : 0;
    }
}

double density_init(int idx) {
    return 1.0;
}

double temperature_init(int idx) {
    auto [i, j, k] = cell_pos(idx);

    if ((i - Nx / 2) * (i - Nx / 2) + (j - Ny / 3) * (j - Ny / 3) < (R) * (R)) {
        return 1.0;
    } else {
        return 0.0;
    }
}


#endif //LBMPOSTELS_LBM_H
