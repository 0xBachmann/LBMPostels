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


void write_vtk(const std::string_view file_name, const std::string& field_name, std::span<const double> data, int Nx, int Ny, double dx);


#endif //LBMPOSTELS_VTK_H
