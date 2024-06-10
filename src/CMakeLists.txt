

add_executable(lbm main.cpp cuda_utils.h constants.h types.h lbm.h lattice.h lattice.h lattice.h vtk.h lbm_cpu.h)

enable_language(CUDA)

add_executable(lbm_gpu main.cu cuda_utils.h constants.h types.h lbm.h lattice.h lattice.h lattice.h vtk.h lbm_cpu.h lbm_gpu.cuh)
target_compile_options(lbm_gpu PUBLIC "--expt-relaxed-constexpr")