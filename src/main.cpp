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


constexpr std::size_t Nx = 200;
constexpr std::size_t Ny = 100;
constexpr std::size_t Nz = 1;
constexpr std::size_t N = Nx*Ny*Nz;
constexpr std::size_t Q = 19;
constexpr std::size_t Npop = N*Q;

// TODO: switch cix dependent on Q
constexpr std::array<double, Q> cqx{0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0};
constexpr std::array<double, Q> cqy{0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
constexpr std::array<double, Q> cqz{0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
constexpr std::array<double, Q> wq{1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};

constexpr double dt = 1.;
constexpr double dx = dt;

// Speed of sound (in lattice units)
constexpr double cs_2 = 3.; // cs^2 = 1./3. * (dx**2/dt**2)
constexpr double cs_4 = 9.;

constexpr double viscosity = 5e-5;
constexpr double relaxation_time = viscosity * cs_2 + 0.5;
constexpr double relaxation_time_inv = 1.0 / relaxation_time;



using PopulationT = std::array<double, Npop>;
using PopulationCellPtr = double*;

using VelocityT = std::array<double, N*3>;
using DensityT = std::array<double, N>;


double f_eq_q(int q, const double ux, const double uy, const double uz, const double u2, const double density) {
    const double uc = ux * cqx[q] + uy * cqy[q] + uz * cqz[q];
    const double uc2 = uc * uc;
    return wq[q] * density * (1. + uc * cs_2 + 0.5 * uc2 * cs_4 - 0.5 * u2 * cs_2);
}

void initialize_population(PopulationT& f_pop, const VelocityT& velocity, const DensityT& density) {
    for(int i = 0; i < N; ++i) {
        const double ux = velocity[i*3];
        const double uy = velocity[i*3 + 1];
        const double uz = velocity[i*3 + 2];
        const double u2 = ux * ux + uy * uy + uz * uz;
        const double cell_density = density[i];
        for(int q = 0; q < Q; ++q) {
            f_pop[i*Q + q] = f_eq_q(q, ux, uy, uz, u2, cell_density);
        }
    }
}

void initialize_moments(VelocityT& velocity, DensityT& density) {
    for(int i = 0; i < N; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            velocity[i*3 + dim] = 0.;
        }
        density[i] = 1.0;     
    }

    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                for(int dim = 0; dim < 3; ++dim) {
                    velocity[i*Ny*Nz*3 + j*Nz*3 + k*3] = 0.0;
                }
            }
        }
    }
}

void relax_population(PopulationT& f_pop, const VelocityT& velocity, const DensityT& density) {
    for(int i = 0; i < N; ++i) {
        const double ux = velocity[i*3];
        const double uy = velocity[i*3 + 1];
        const double uz = velocity[i*3 + 2];
        const double u2 = ux * ux + uy * uy + uz * uz;
        const double cell_density = density[i];
        for(int q = 0; q < Q; ++q) {
            f_pop[i*Q + q] = (1.0 - relaxation_time_inv) * f_pop[i*Q + q] + relaxation_time_inv * f_eq_q(q, ux, uy, uz, u2, cell_density);
        }
    }
}

int cell_pop_idx(int i, int j, int k) {
    return (i%Nx)*Ny*Nz*Q + (j%Ny)*Nz*Q + (k%Nz)*Q;
}

// BB Boundary
// If q is odd we have to add one to it to find the flipped vector index
// If even we remove one
template <typename P>
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
    std::swap(f_pop, f_pop_buffer);
    // Boundaries
    // for(int j = 1; j < Ny - 1 ; ++j){
    //     for(int k = 1; k < Nz - 1; ++k) {
    //         // Left x-
    //         int i = 0;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqx[q] == 1;});

    //         // Right x+
    //         i = Nx - 1;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqx[q] == -1;});
    //     }
    // }
    // for(int i = 1; i < Nx - 1; ++i) {
    //     for(int k = 1; k < Nz - 1; ++k) {
    //         // Bottom: y-
    //         int j = 0;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqy[q] == 1;});

    //         // Top: y+
    //         j = Ny - 1;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqy[q] == -1;});
    //     }
    
    //     for(int j = 1; j < Ny; ++j) {
    //         // Back: z-
    //         int k = 0;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqz[q] == 1;});

    //         // Front: z+
    //         k = Nz - 1;
    //         apply_BB_boundary(f_pop, f_pop_buffer, i, j, k, [](int q){return cqz[q] == -1;});
    //     }
    // }

    // // Corners
    // // Left Bottom Back
    // apply_BB_boundary(f_pop, f_pop_buffer, 0, 0, 0, [](int q){return cqx[q] == 1 || ;});


}

void update_moments(VelocityT& velocity, DensityT& density, const PopulationT& f_pop) {
    for(int i = 0; i < N; ++i) {
        double cell_density = 0.;
        double ux = 0.;
        double uy = 0.;
        double uz = 0.;

        for(int q = 0; q < Q; ++q) {
            cell_density += f_pop[i*Q + q];
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
    }
}

void apply_force(VelocityT& velocity) {
    // Gravity?
    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = 0; k < Nz; ++k) {
                velocity[i*Ny*Nz*3 + j*Nz*3 + k*3] += 0.001;
            }
        }
    }
    for(int i = 2*Nx/5; i < 3*Nx/5; ++i){
        for(int j = 2*Ny/5; j < 3*Ny/5; ++j) {
            for(int k = 0; k < Nz; ++k) {
                for(int dim = 0; dim < 3; ++dim) {
                    velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim] = 0.99*velocity[i*Ny*Nz*3 + j*Nz*3 + k*3 + dim];
                }
            }
        }
    }
}
void write_vtk(const std::string_view file_name, const DensityT& density, const VelocityT& velocity) {
    std::ofstream file{std::string(file_name)};

    file << "# vtk DataFile Version 3.0\n";
    file << "LBM DATA" << "\n"; // Title
    file << "ASCII\n";
    
    file << std::format("DATASET STRUCTURED_POINTS\n"
                        "DIMENSIONS {} {} {}\n"
                        "ORIGIN 0 0 0\n"
                        "SPACING {} {} {}\n\n"
                        "POINT_DATA {}\n", Nx, Ny, Nz, dx, dx, dx, N);
    // constexpr std::string scalars_data_string = "SCALARS {} double 1\n"
    //                                   "LOOKUP_TABLE default\n";
    
    file << "VECTORS velocity double\n";
    for (int i = 0; i < N; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            file << velocity[i*3 + dim] << " ";
        }
        file << "\n";
    }

    // file << std::format(scalars_data_string, "velocity_y");
    // for (int i = 0; i < N; ++i) {
    //     file << velocity[i*3 + 1] << " ";
    // }
    // file << "\n";

    // file << std::format(scalars_data_string, "velocity_z");
    // for (int i = 0; i < N; ++i) {
    //     file << velocity[i*3 + 2] << " ";
    // }
    // file << "\n";

    // file << std::format(scalars_data_string, "density");
    // for (int i = 0; i < N; ++i) {
    //     file << density[i] << " ";
    // }
    // file << "\n";
}

int main() {
    PopulationT f_pop;
    PopulationT f_pop_buffer;

    VelocityT velocity;
    DensityT density;

    // TODO: make this user input
    const int Nsteps = 1000;

    // 1. Initial conditions 
    initialize_moments(velocity, density);
    initialize_population(f_pop, velocity, density);

    // 2. Time stepping
    for(int t = 0; t < Nsteps; ++t){
        if(!((t+1) % 10))
            std::cout << std::format("{}/{} steps \n", t+1, Nsteps);
        // 2.0 Forcing
        apply_force(velocity);
        // 2.1 Relaxation
        relax_population(f_pop, velocity, density);
        // 2.2 Streaming
        stream_population(f_pop, f_pop_buffer);
        // 2.3 Update moments
        update_moments(velocity, density, f_pop);
        if(!(t % 10))
            write_vtk(std::format("hello-{}.vtk",t/10) , density, velocity);

    }


    write_vtk("hello.vtk", density, velocity);
}