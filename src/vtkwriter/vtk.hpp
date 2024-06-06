/*!
 * @author jonas 
 * @date 05.06.24
 */

#ifndef LBMPOSTELS_VTK_H
#define LBMPOSTELS_VTK_H

#include <format>
#include <string_view>
#include <iostream>
#include <fstream>
#include <ostream>

#include <Eigen/Eigen>

void write_vtk(std::string_view file_name, const Eigen::Matrix<std::array<double, 2>, -1, -1>& grid, double dx);


#endif //LBMPOSTELS_VTK_H
