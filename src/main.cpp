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

#include "lbm_cpu.h"

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
    for (int t = 0; t < Nsteps; ++t) {
        if (!((t + 1) % 10))
            std::cout << std::format("{}/{} steps\r", t + 1, Nsteps) << std::endl;
        // 2.0 Forcing
        apply_force(velocity);

        // 2.1 Relaxation
        relax_populations(f_pop, g_pop, velocity, density, temperature);
        // 2.2 Streaming
        stream_population(f_pop, f_pop_buffer);
        stream_population(g_pop, g_pop_buffer);

        // 2.3 Update moments
        update_moments(velocity, density, temperature, f_pop, g_pop);
        if (!(t % output_frequency))
            write_vtk(std::format("hello-{}.vtk", t / output_frequency + 1), density, velocity, temperature);
    }

    write_vtk("hello.vtk", density, velocity, temperature);
}