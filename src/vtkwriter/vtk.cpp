/*!
 * @author jonas 
 * @date 05.06.24
 */

#include "vtk.hpp"

void write_vtk(const std::string_view file_name, const std::string& field_name, std::span<const double> data, int Nx, int Ny, double dx) {
    std::ofstream file{std::string(file_name)};

    file << "# vtk DataFile Version 3.0\n";
    file << field_name << "\n"; // Title
    file << "ASCII\n";
    file << std::format("DATASET STRUCTURED_POINTS\n"
                        "DIMENSIONS {} {} 1\n"
                        "ORIGIN 0 0 0\n"
                        "SPACING {} {} 1\n\n"
                        "POINT_DATA {}\n"
                        "SCALARS {} double 1\n"
                        "LOOKUP_TABLE default\n", Nx, Ny, dx, dx, Nx*Ny, field_name);
                        
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            // file << std::format("{} ", data[i*Ny + j]);
            file << 10.0 << " ";

        }
    }
}
