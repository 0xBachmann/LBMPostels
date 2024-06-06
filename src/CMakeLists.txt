find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_executable(lbm vtkwriter/vtk.cpp main.cpp vtkwriter/vtk.hpp grid.hpp)
target_link_libraries(lbm Eigen3::Eigen)