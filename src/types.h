#ifndef LBMPOSTELS_TYPES_H
#define LBMPOSTELS_TYPES_H

#include <array>

#include "constants.h"

#ifdef __NVCC__
#include <thrust/universal_vector.h>

using PopulationT = thrust::universal_vector<double>;
using PopulationCellPtr = double*;

using VelocityT = thrust::universal_vector<double>;
using ScalarT = thrust::universal_vector<double>;
#else
using PopulationT = std::array<double, Npop>;
using PopulationCellPtr = double*;

using VelocityT = std::array<double, N*3>;
using ScalarT = std::array<double, N>;
#endif

#endif //LBMPOSTELS_TYPES_H
