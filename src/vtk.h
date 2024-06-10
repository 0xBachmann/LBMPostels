#ifndef LBMPOSTELS_VTK_H
#define LBMPOSTELS_VTK_H

#include <string>
#include <string_view>
#include <format>
#include <iostream>
#include <fstream>

#include "types.h"
#include "constants.h"

void write_vtk(const std::string_view file_name, std::span<const double> density, std::span<const double> velocity, std::span<const double> temperature) {
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


#endif //LBMPOSTELS_VTK_H
