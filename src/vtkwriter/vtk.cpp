/*!
 * @author jonas 
 * @date 05.06.24
 */

#include "vtk.hpp"

void write_vtk(std::string_view file_name, const Eigen::Matrix<std::array<double, 2>, -1, -1>& grid, double dx) {
    std::ofstream file{std::string(file_name)};

    file << "# vtk DataFile Version 2.0\n";
    file << "\n"; // comment
    file << "ASCII\n";
    file << std::format("DATASET STRUCTURED_POINTS\n"
                        "DIMENSIONS {} {} 0\n"
                        "ORIGIN 0 0 0\n"
                        "SPACING {} {} 1\n"
                        "CELL_DATA {}\n"
                        "VECTORS velocity double\n", grid.rows(), grid.cols(), dx, dx, grid.rows() * grid.cols());
    for (const auto& col : grid.colwise()) {
        for (const auto& vel : col) {
            file << std::format("{}, {}\n", vel[0], vel[1]);
        }
    }
}
