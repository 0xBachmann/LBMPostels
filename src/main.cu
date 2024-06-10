#include <iostream>

#include "lbm_gpu.cuh"
#include "vtk.h"

int main() {
    thrust::universal_vector<double> f_pop_v(Npop);
    thrust::universal_vector<double> f_pop_buffer_v(Npop);

    thrust::universal_vector<double> g_pop_v(Npop);
    thrust::universal_vector<double> g_pop_buffer_v(Npop);

    thrust::universal_vector<double> velocity_v(N * 3);
    thrust::universal_vector<double> density_v(N);
    thrust::universal_vector<double> temperature_v(N);

    std::span<double> f_pop{thrust::raw_pointer_cast(f_pop_v.data()), f_pop_v.size()};
    std::span<double> f_pop_buffer{thrust::raw_pointer_cast(f_pop_buffer_v.data()), f_pop_buffer_v.size()};

    std::span<double> g_pop{thrust::raw_pointer_cast(g_pop_v.data()), g_pop_v.size()};
    std::span<double> g_pop_buffer{thrust::raw_pointer_cast(g_pop_buffer_v.data()), g_pop_buffer_v.size()};

    std::span<double> velocity{thrust::raw_pointer_cast(velocity_v.data()), velocity_v.size()};
    std::span<double> density{thrust::raw_pointer_cast(density_v.data()), density_v.size()};
    std::span<double> temperature{thrust::raw_pointer_cast(temperature_v.data()), temperature_v.size()};

    // TODO: make this user input
    const int Nsteps = 10000;

    constexpr int num_threads = 512;
    constexpr int num_blocks_n = (N + num_threads - 1) / num_threads;
    constexpr int num_blocks_nq = (N*Q + num_threads - 1) / num_threads;

    // 1. Initial conditions
    initialize_moments_gpu<<<num_blocks_n, num_threads>>>(velocity, density, temperature);
    initialize_populations_gpu<<<num_blocks_n, num_threads>>>(f_pop, g_pop, velocity, density, temperature);

    // 2. Time stepping
    for (int t = 0; t < Nsteps; ++t) {
        if (!((t + 1) % 10))
            std::cout << std::format("{}/{} steps\r", t + 1, Nsteps) << std::endl;
        // 2.0 Forcing
        apply_force_gpu<<<num_blocks_n, num_threads>>>(velocity);

        // 2.1 Relaxation
        relax_populations_gpu<<<num_blocks_nq, num_threads>>>(f_pop, g_pop, velocity, density, temperature);
        // 2.2 Streaming
        stream_population_gpu<<<num_blocks_nq, num_threads>>>(f_pop, f_pop_buffer);
        f_pop_v = f_pop_buffer_v;
        stream_population_gpu<<<num_blocks_nq, num_threads>>>(g_pop, g_pop_buffer);
        g_pop_v = g_pop_buffer_v;

        // 2.3 Update moments
        update_moments_gpu<<<num_blocks_n, num_threads>>>(velocity, density, temperature, f_pop, g_pop);
        if (!(t % output_frequency))
            write_vtk(std::format("hello-{}.vtk", t / output_frequency + 1), density, velocity, temperature);
    }

    write_vtk("hello.vtk", density, velocity, temperature);
}