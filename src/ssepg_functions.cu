/*
  This file is part of the MRF_CUDA package (https://github.com/chixindebaoyu/MRF_CUDA)

  The MIT License (MIT)

  Copyright (c) 2019 Dong Wang and David Smith

  Permission is hereby granted, free of charge, to any person obtaining a # copy
  of this software and associated documentation files (the "Software"), to # deal
  in the Software without restriction, including without limitation the # rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or # sell
  coM_PIes of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included # in all
  coM_PIes or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS # OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL # THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING # FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS # IN THE
  SOFTWARE.
*/

// #include "functions.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "cuda_runtime.h"
#include <cublas_v2.h>
#include "float2math.h"

// Specify the sizes according to the divice
// int blocksize = 64;
// int gridsize = 128;

// kernels


__global__ void
dephase_gradients_rf_stage1(float2 * d_w_out, float2 * d_w, size_t natoms, int nstates)
{
    int ii, jj, idx;
    int from_idx, to_idx;

    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < 2 * natoms * nstates - 2;
        id += blockDim.x * gridDim.x)
    {

		ii = id / 2;
		jj = natoms * nstates - ii - 1;
        idx = (id + 1) % 2;

        from_idx = (1 - idx) * jj + idx * ii;
        to_idx = from_idx + 2 * idx - 1;

        d_w_out[3 * from_idx + 1 - idx] =
            d_w[3 * to_idx + 1 - idx];
    }
}

__global__ void
dephase_gradients_rf_stage2(float2 * d_w, size_t natoms, int nstates)
{
    int idx;

    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < natoms;
        id += blockDim.x * gridDim.x)
    {
        idx = id * 3 * nstates + 1;

        d_w[idx] = conj(d_w[idx - 1]);
        d_w[idx + 3 * nstates - 4] = make_float2(0.f);
    }
}

