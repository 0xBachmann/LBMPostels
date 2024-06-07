/*!
 * @author jonas 
 * @date 05.06.24
 */

#include <array>
#include <format>
#include <string_view>
#include <iostream>
#include <fstream>
#include <ostream>

#include "constants.h"
#include "types.h"
#include "lbm.h"

void initialize_populations(PopulationT& f_pop, PopulationT& g_pop, const VelocityT& velocity, const ScalarT& density, const ScalarT& temperature) {
    for(int i = 0; i < N; ++i) {
        auto init_data = get_cell_init_data(i, velocity, density, temperature);
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
            //velocity[i*3 + dim] = 0.;
        }
        //velocity[i*3 + 1] = U;
        density[i] = 1.0;
        temperature[i] = temperature_init(i);
    }

    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                // velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + 1] += 0.001;
                //velocity[i*Ny*Nz*3 + j*Nz*3 + k*3+1] += 0.001;
                if((i - Nx/2)*(i - Nx/2) + (j - Ny/3)*(j - Ny/3) < (R)*(R)) {
                    for(int dim = 0; dim < 3; ++dim) {
                        // velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] = 0;
                        //velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] = 0;
                        // temperature[i*Ny*Nz + j*Nz + k] = 1.0;
                    }
                }
            }
        }
    }
}

void relax_populations(PopulationT& f_pop, PopulationT& g_pop, const VelocityT& velocity, const ScalarT& density, const ScalarT& temperature) {
    for(int i = 0; i < N; ++i) {
        const double ux = velocity[i*3];
        const double uy = velocity[i*3 + 1];
        const double uz = velocity[i*3 + 2];
        const double u2 = ux * ux + uy * uy + uz * uz;
        const double cell_density = density[i];
        const double cell_temperature = temperature[i];
        for(int q = 0; q < Q; ++q) {
            f_pop[i*Q + q] = (1.0 - relaxation_time_inv) * f_pop[i*Q + q] + relaxation_time_inv * f_eq_q(q, ux, uy, uz, u2, cell_density);
            g_pop[i*Q + q] = (1.0 - relaxation_time_g_inv) * g_pop[i*Q + q] + relaxation_time_g_inv * g_eq_q(q, ux, uy, uz, u2, cell_temperature);
        }
    }
}

// BB Boundary
// If q is odd we have to add one to it to find the flipped vector index
// If even we remove one
template<typename P>
void apply_BB_boundary(PopulationT& f_pop, PopulationT& f_pop_buffer, int i, int j, int k, P&& predicate) {
    const double current_cell_idx = cell_pop_idx(i, j, k);
    for(int q = 0; q < Q; ++q) {
        if(predicate(q)) {
            f_pop_buffer[current_cell_idx + q] = f_pop[current_cell_idx + (q%2 ? q + 1 : q - 1)];
        } else {
            f_pop_buffer[current_cell_idx + q] = f_pop[cell_pop_idx(i - cqx[q], j - cqy[q], k - cqz[q]) + q];
        }
    }
}

void stream_population(PopulationT& f_pop, PopulationT& f_pop_buffer) {
    // Inner cells
    for(int i = 0; i < N; ++i) {
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
        double cell_density = 0.;
        double ux = 0.;
        double uy = 0.;
        double uz = 0.;
        double cell_temperature = 0.;

        for(int q = 0; q < Q; ++q) {
            cell_density += f_pop[i*Q + q];
            cell_temperature += g_pop[i*Q + q];
            ux += cqx[q] * f_pop[i*Q + q];
            uy += cqy[q] * f_pop[i*Q + q];
            uz += cqz[q] * f_pop[i*Q + q]; 
        }

        double cell_density_inv = 1./cell_density;
        ux *= cell_density_inv;
        uy *= cell_density_inv;
        uz *= cell_density_inv;

        velocity[i*3] = ux;
        velocity[i*3 + 1] = uy;
        velocity[i*3 + 2] = uz;

        density[i] = cell_density;
        temperature[i] = cell_temperature;
    }
}

void apply_force(VelocityT& velocity) {
    // Gravity?
    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                // velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + 1] += 0.001;
                //velocity[i*Ny*Nz*3 + j*Nz*3 + k*3+1] += 0.001;
                // if(j < Ny/5)
                    // velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + 1] += 0.001;
                if((i - Nx/2)*(i - Nx/2) + (j - Ny/3)*(j - Ny/3) < (R)*(R)) {
                    for(int dim = 0; dim < 3; ++dim) {
                        // velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] = 0.9*velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim];
                        velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] = 0;
                    }                    
                }
            }
        }
    }
}
void write_vtk(const std::string_view file_name, const ScalarT& density, const VelocityT& velocity, const ScalarT& temperature) {
    std::ofstream file{std::string(file_name)};

    file << "# vtk DataFile Version 3.0\n";
    file << "LBM DATA\n"; // Title
    file << "ASCII\n";
    
    file << std::format("DATASET STRUCTURED_POINTS\n"
                        "DIMENSIONS {} {} {}\n"
                        "ORIGIN 0 0 0\n"
                        "SPACING {} {} {}\n\n"
                        "POINT_DATA {}\n", Nx, Ny, Nz, dx, dx, dx, N);
    
    file << "VECTORS velocity double\n";
    for(int k = 0; k < Nz; ++k) {
        for(int j = 0; j < Ny; ++j) {
            for(int i = 0; i < Nx; ++i) {
                for(int dim = 0; dim < 3; ++dim) {
                    file << velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] << " ";
                }
                file << "\n";
            }
        }
    }

    file << "\n";
    file << "SCALARS density double 1\n"
            "LOOKUP_TABLE default\n";
    for(int k = 0; k < Nz; ++k) {
        for(int j = 0; j < Ny; ++j) {
            for(int i = 0; i < Nx; ++i) {
                file << density[i*Ny*Nz + j*Nz + k] << " ";
            }
        }
    }
    file << "\n\n";

    file << "SCALARS temperature double 1\n"
            "LOOKUP_TABLE default\n";
    for(int k = 0; k < Nz; ++k) {
        for(int j = 0; j < Ny; ++j) {
            for(int i = 0; i < Nx; ++i) {
                file << temperature[i*Ny*Nz + j*Nz + k] << " ";
            }
        }
    }
    file << "\n";
}

int main() {
    PopulationT f_pop;
    PopulationT f_pop_buffer;

    PopulationT g_pop;
    PopulationT g_pop_buffer;

    VelocityT velocity;
    ScalarT density;
    ScalarT temperature;

    // TODO: make this user input
    const int Nsteps = 10000;

    // 1. Initial conditions 
    initialize_moments(velocity, density, temperature);
    initialize_populations(f_pop, g_pop, velocity, density, temperature);

    // 2. Time stepping
    for(int t = 0; t < Nsteps; ++t){
        if(!((t+1) % 10))
            std::cout << std::format("{}/{} steps\r", t+1, Nsteps) << std::endl;
        // 2.0 Forcing
        apply_force(velocity);

        // 2.1 Relaxation
        relax_populations(f_pop, g_pop, velocity, density, temperature);
        // 2.2 Streaming
        stream_population(f_pop, f_pop_buffer);
        stream_population(g_pop, g_pop_buffer);

        // 2.3 Update moments
        update_moments(velocity, density, temperature, f_pop, g_pop);
        if(!(t % output_frequency))
            write_vtk(std::format("hello-{}.vtk",t/output_frequency+1) , density, velocity, temperature);
    }


    write_vtk("hello.vtk", density, velocity, temperature);
}