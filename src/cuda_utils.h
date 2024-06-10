#ifndef LBMPOSTELS_CUDA_UTILS_H
#define LBMPOSTELS_CUDA_UTILS_H

#ifdef __NVCC__
#define CUDA_HOST_DEVICE __host__ __device__
#else
#define CUDA_HOST_DEVICE
#endif

#ifdef __NVCC__
#define CUDA_CONST __constant__
#else
#define CUDA_CONST
#endif

#endif //LBMPOSTELS_CUDA_UTILS_H
