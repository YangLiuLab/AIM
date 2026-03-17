/*=================================================================
 * PointIntersect1D - Point-based 1D cross-correlation for Z drift (OpenMP)
 *
 * Pure C++ MEX source file. Compile with:
 *   Windows (MSVC):  mex PointIntersect1D.cpp COMPFLAGS="$COMPFLAGS /openmp"
 *   macOS (Xcode clang default has no OpenMP, use Homebrew gcc):
 *       mex PointIntersect1D.cpp CXX='/opt/homebrew/bin/g++-14' CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *   Linux (GCC):     mex PointIntersect1D.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
 *
 * If OpenMP is unavailable, compile without flags:
 *       mex PointIntersect1D.cpp
 *
 * Usage in MATLAB:
 *   ROIcc = PointIntersect1D(Pxyz0, Vxyz0, length(Vxyz0), Pxyz1, Vxyz1, length(Vxyz1), imW/IntersectD, roiR)
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

static void cross_corr_1D_omp(
    double* plist_1, double* pval_1, int pnum_1,
    double* plist_n, double* pval_n, int pnum_n,
    double* img_crr, double step, int roiR)
{
    int sz_corr = roiR * 2 + 1;
    long long z_stride = (long long)(step * step);

#ifdef _OPENMP
    const int g_ncore = omp_get_num_procs();
    int tn = sz_corr > g_ncore ? g_ncore : sz_corr;
    if (tn < 1) tn = 1;
#else
    int tn = 1;
#endif

#pragma omp parallel for schedule(dynamic) num_threads(tn)
    for (int ri = 0; ri < sz_corr; ri++)
    {
        int r = ri - roiR;  /* shift from -roiR to +roiR */
        long long sft = (long long)r * z_stride;
        int temp = 0;
        int sidx = 0;

        for (int pidx_n = 0; pidx_n < pnum_n; pidx_n++)
        {
            long long cur_n = (long long)plist_n[pidx_n] + sft;

            for (int pidx_1 = sidx; pidx_1 < pnum_1; pidx_1++)
            {
                if (cur_n == (long long)plist_1[pidx_1])
                {
                    temp += (int)(pval_1[pidx_1] * pval_n[pidx_n]);
                    sidx = pidx_1 + 1;
                    break;
                }
                else if (cur_n < (long long)plist_1[pidx_1])
                {
                    break;
                }
            }
        }
        img_crr[ri] = (double)temp;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 8)
        mexErrMsgIdAndTxt("PointIntersect1D:nrhs", "Eight inputs required.");

    double* plist_1 = mxGetPr(prhs[0]);
    double* pval_1  = mxGetPr(prhs[1]);
    int     pnum_1  = (int)mxGetScalar(prhs[2]);
    double* plist_n = mxGetPr(prhs[3]);
    double* pval_n  = mxGetPr(prhs[4]);
    int     pnum_n  = (int)mxGetScalar(prhs[5]);
    double  step    = mxGetScalar(prhs[6]);
    int     roiR    = (int)mxGetScalar(prhs[7]);

    int sz_corr = roiR * 2 + 1;

    plhs[0] = mxCreateDoubleMatrix(sz_corr, 1, mxREAL);
    double* img_crr = mxGetPr(plhs[0]);

    cross_corr_1D_omp(plist_1, pval_1, pnum_1,
                       plist_n, pval_n, pnum_n,
                       img_crr, step, roiR);
}
