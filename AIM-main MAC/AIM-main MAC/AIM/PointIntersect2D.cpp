/*=================================================================
 * PointIntersect2D - Point-based 2D cross-correlation (OpenMP)
 *
 * Pure C++ MEX source file. Compile with:
 *   Windows (MSVC):  mex PointIntersect2D.cpp COMPFLAGS="$COMPFLAGS /openmp"
 *   macOS (Xcode clang default has no OpenMP, use Homebrew gcc):
 *       mex PointIntersect2D.cpp CXX='/opt/homebrew/bin/g++-14' CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *   Linux (GCC):     mex PointIntersect2D.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
 *
 * If OpenMP is unavailable, compile without flags:
 *       mex PointIntersect2D.cpp
 *
 * Usage in MATLAB:
 *   img_crr = PointIntersect2D(Pxy0, Vxy0, length(Vxy0), Pxy1, Vxy1, length(Vxy1), array_sft, array_idx, array_len, roiR)
 *
 *=================================================================*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include "mex.h"
#include "matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static void cross_corr_2D_omp(
    double* plist_1, double* pval_1, int pnum_1,
    double* plist_n, double* pval_n, int pnum_n,
    double* img_crr, double* array_sft, double* array_idx, int array_len)
{
    int i = 0;

#ifdef _OPENMP
    const int g_ncore = omp_get_num_procs();
    int max_tn = array_len / 4;
    int tn = max_tn > g_ncore ? g_ncore : max_tn;
    if (tn < 1) tn = 1;
#else
    int tn = 1;
#endif

#pragma omp parallel for schedule(dynamic) num_threads(tn)
    for (i = 0; i < array_len; i++)
    {
        int sft = (int)array_sft[i];
        int temp = 0;
        int sidx = 0;

        for (int pidx_n = 0; pidx_n < pnum_n; pidx_n++)
        {
            int cur_n = (int)plist_n[pidx_n] + sft;

            for (int pidx_1 = sidx; pidx_1 < pnum_1; pidx_1++)
            {
                if (cur_n == (int)plist_1[pidx_1])
                {
                    temp += (int)(pval_1[pidx_1] * pval_n[pidx_n]);
                    sidx = pidx_1 + 1;
                    break;
                }
                else if (cur_n < (int)plist_1[pidx_1])
                {
                    break;
                }
            }
        }
        img_crr[(int)array_idx[i]] = (double)temp;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 10)
        mexErrMsgIdAndTxt("PointIntersect2D:nrhs", "Ten inputs required.");

    double* plist_1   = mxGetPr(prhs[0]);
    double* pval_1    = mxGetPr(prhs[1]);
    int     pnum_1    = (int)mxGetScalar(prhs[2]);
    double* plist_n   = mxGetPr(prhs[3]);
    double* pval_n    = mxGetPr(prhs[4]);
    int     pnum_n    = (int)mxGetScalar(prhs[5]);
    double* array_sft = mxGetPr(prhs[6]);
    double* array_idx = mxGetPr(prhs[7]);
    int     array_len = (int)mxGetScalar(prhs[8]);
    int     r_corr    = (int)mxGetScalar(prhs[9]);

    int sz_corr  = r_corr * 2 + 1;
    int len_corr = sz_corr * sz_corr;

    plhs[0] = mxCreateDoubleMatrix(len_corr, 1, mxREAL);
    double* img_crr = mxGetPr(plhs[0]);

    cross_corr_2D_omp(plist_1, pval_1, pnum_1,
                       plist_n, pval_n, pnum_n,
                       img_crr, array_sft, array_idx, array_len);
}
