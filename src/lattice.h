#ifndef LBMPOSTELS_LATTICE_H
#define LBMPOSTELS_LATTICE_H

#include <array>

#include "constants.h"

int cell_pop_idx(int i, int j, int k) {
    return ((i+Nx)%Nx)*Ny*Nz*Q + ((j+Ny)%Ny)*Nz*Q + ((k+Nz)%Nz)*Q;
}


// returns (i,j,k)
std::array<int, 3> cell_pos(int idx) {
    int k = idx % Nz;
    int j = (idx / Nz) % Ny;
    int i = (idx / (Nz * Ny));
    return {i,j,k};
}


// returns (i,j,k,q)
std::array<int, 4> cell_vec_pos(int idx) {
    int q = idx % Q;
    int k = (idx / Q) % Nz;
    int j = (idx / (Q * Nz)) % Ny;
    int i = (idx / (Q * Nz * Ny));
    return {i,j,k,q};
}

#endif //LBMPOSTELS_LATTICE_H
