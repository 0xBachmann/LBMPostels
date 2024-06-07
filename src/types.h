#ifndef LBMPOSTELS_TYPES_H
#define LBMPOSTELS_TYPES_H

#include <array>

#include "constants.h"

using PopulationT = std::array<double, Npop>;
using PopulationCellPtr = double*;

using VelocityT = std::array<double, N*3>;
using ScalarT = std::array<double, N>;

#endif //LBMPOSTELS_TYPES_H
