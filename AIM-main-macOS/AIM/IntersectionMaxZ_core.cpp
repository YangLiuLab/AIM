/*=================================================================
 * IntersectionMaxZ_core - Optimized Z drift tracking loop
 *
 * Same optimizations as IntersectionMax_core:
 *   1. Pre-bin points by segment upfront
 *   2. Binary search in cross-correlation
 *   3. OpenMP parallel shifts
 *
 * Usage:
 *   driftZ = IntersectionMaxZ_core(XList, YList, ZList, refXList, refYList,
 *               refZList, fID, trackNUM, trackInterval, imW, IntersectD)
 *=================================================================*/
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "mex.h"
#include "matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void groupcounts_sorted(
    const std::vector<long long>& sorted,
    std::vector<long long>& values,
    std::vector<int>& counts)
{
    values.clear();
    counts.clear();
    if (sorted.empty()) return;
    values.reserve(sorted.size() / 2);
    counts.reserve(sorted.size() / 2);
    values.push_back(sorted[0]);
    counts.push_back(1);
    for (size_t i = 1; i < sorted.size(); i++) {
        if (sorted[i] == values.back()) {
            counts.back()++;
        } else {
            values.push_back(sorted[i]);
            counts.push_back(1);
        }
    }
}

static inline int binary_find(const long long* arr, int n, long long key) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (arr[mid] == key) return mid;
        else if (arr[mid] < key) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}

static void cross_corr_1d_fast(
    const long long* ref_pos, const int* ref_val, int ref_n,
    const long long* tgt_pos, const int* tgt_val, int tgt_n,
    long long z_stride, int roiR, double* output)
{
    int sz = 2 * roiR + 1;

#ifdef _OPENMP
    int ncores = omp_get_num_procs();
    int tn = sz < ncores ? sz : ncores;
    if (tn < 1) tn = 1;
#endif

    #pragma omp parallel for schedule(static) num_threads(tn)
    for (int ri = 0; ri < sz; ri++) {
        long long sft = (long long)(ri - roiR) * z_stride;
        long long temp = 0;

        for (int j = 0; j < tgt_n; j++) {
            long long key = tgt_pos[j] + sft;
            int idx = binary_find(ref_pos, ref_n, key);
            if (idx >= 0) {
                temp += (long long)ref_val[idx] * tgt_val[j];
            }
        }
        output[ri] = (double)temp;
    }
}

static double dft_peak_1d(const double* data, int N) {
    double freq = 2.0 * M_PI / (double)N;
    double re = 0.0, im = 0.0;
    for (int k = 0; k < N; k++) {
        re += data[k] * cos(freq * k);
        im -= data[k] * sin(freq * k);
    }
    double ang = atan2(im, re);
    if (ang > 0) ang -= 2.0 * M_PI;
    return (fabs(ang) / (2.0 * M_PI / N) + 1.0) - (N + 1.0) / 2.0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 11)
        mexErrMsgIdAndTxt("IntersectionMaxZCore:nrhs", "Eleven inputs required.");

    double* XList    = mxGetPr(prhs[0]);
    double* YList    = mxGetPr(prhs[1]);
    double* ZList    = mxGetPr(prhs[2]);
    double* refXList = mxGetPr(prhs[3]);
    double* refYList = mxGetPr(prhs[4]);
    double* refZList = mxGetPr(prhs[5]);
    double* fID_d    = mxGetPr(prhs[6]);
    int trackNUM     = (int)mxGetScalar(prhs[7]);
    int trackInterval = (int)mxGetScalar(prhs[8]);
    double imW       = mxGetScalar(prhs[9]);
    double IntersectD = mxGetScalar(prhs[10]);

    int nPoints = (int)mxGetNumberOfElements(prhs[0]);
    int nRef    = (int)mxGetNumberOfElements(prhs[3]);

    double scale = imW / IntersectD;
    long long z_stride = (long long)(scale * scale);
    int roiR = 3;
    int ROI_size = 2 * roiR + 1;

    /* Pre-encode all positions */
    std::vector<long long> pList(nPoints);
    std::vector<int> fID(nPoints);
    for (int i = 0; i < nPoints; i++) {
        pList[i] = (long long)round(ZList[i] / IntersectD) * z_stride
                 + (long long)round(YList[i] / IntersectD) * (long long)scale
                 + (long long)round(XList[i] / IntersectD);
        fID[i] = (int)fID_d[i];
    }

    /* Pre-bin by segment using sorted index */
    std::vector<int> sortIdx(nPoints);
    for (int i = 0; i < nPoints; i++) sortIdx[i] = i;
    std::sort(sortIdx.begin(), sortIdx.end(),
        [&fID](int a, int b) { return fID[a] < fID[b]; });

    std::vector<int> seg_start(trackNUM + 1, nPoints);
    std::vector<int> seg_end(trackNUM + 1, 0);
    for (int i = 0; i < nPoints; i++) {
        int f = fID[sortIdx[i]];
        int seg = (f - 1) / trackInterval;
        if (seg >= trackNUM) seg = trackNUM - 1;
        if (i < seg_start[seg]) seg_start[seg] = i;
        if (i + 1 > seg_end[seg]) seg_end[seg] = i + 1;
    }

    /* Encode reference and groupcounts */
    std::vector<long long> refEncoded(nRef);
    for (int i = 0; i < nRef; i++) {
        refEncoded[i] = (long long)round(refZList[i] / IntersectD) * z_stride
                      + (long long)round(refYList[i] / IntersectD) * (long long)scale
                      + (long long)round(refXList[i] / IntersectD);
    }
    std::sort(refEncoded.begin(), refEncoded.end());

    std::vector<long long> ref_vals;
    std::vector<int> ref_cnts;
    groupcounts_sorted(refEncoded, ref_vals, ref_cnts);
    int ref_n = (int)ref_vals.size();

    /* Output */
    plhs[0] = mxCreateDoubleMatrix(1, trackNUM, mxREAL);
    double* driftZ = mxGetPr(plhs[0]);

    /* Pre-allocate buffers */
    std::vector<long long> seg_encoded;
    std::vector<long long> seg_vals;
    std::vector<int> seg_cnts;
    std::vector<double> ROIcc(ROI_size);
    seg_encoded.reserve(nPoints / trackNUM * 2);

    /* Main loop */
    double refz = 0.0;

    for (int s = 0; s < trackNUM; s++) {
        seg_encoded.clear();
        if (seg_start[s] < seg_end[s]) {
            for (int i = seg_start[s]; i < seg_end[s]; i++) {
                int idx = sortIdx[i];
                if (fID[idx] > s * trackInterval && fID[idx] <= (s + 1) * trackInterval) {
                    seg_encoded.push_back(pList[idx]);
                }
            }
        }

        if (seg_encoded.empty()) {
            driftZ[s] = (s > 0) ? driftZ[s-1] : 0.0;
            continue;
        }

        std::sort(seg_encoded.begin(), seg_encoded.end());
        groupcounts_sorted(seg_encoded, seg_vals, seg_cnts);
        int seg_n = (int)seg_vals.size();

        long long sft = (long long)round(refz) * z_stride;
        for (int i = 0; i < seg_n; i++) {
            seg_vals[i] += sft;
        }

        std::fill(ROIcc.begin(), ROIcc.end(), 0.0);
        cross_corr_1d_fast(
            ref_vals.data(), ref_cnts.data(), ref_n,
            seg_vals.data(), seg_cnts.data(), seg_n,
            z_stride, roiR, ROIcc.data());

        double PZ = dft_peak_1d(ROIcc.data(), ROI_size);

        refz = round(refz) + PZ;
        driftZ[s] = -refz;
    }
}