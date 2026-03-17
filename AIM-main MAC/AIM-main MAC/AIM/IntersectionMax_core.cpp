/*=================================================================
 * IntersectionMax_core - Optimized full drift tracking loop
 *
 * Key optimizations vs previous version:
 *   1. Pre-bins all points by segment ONCE upfront (O(N) total vs O(N*trackNUM))
 *   2. Binary search in cross-correlation (O(log N) vs O(N) per lookup)
 *   3. Pre-sorted reference list, no re-sorting per iteration
 *   4. Zero allocations inside the hot loop
 *   5. OpenMP parallel cross-correlation
 *
 * Compile: see compile_mex.m
 *
 * Usage:
 *   [driftX, driftY] = IntersectionMax_core(XList, YList, refXList, refYList,
 *                          fID, trackNUM, trackInterval, imW, IntersectD)
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

/*---------------------------------------------------------------
 * Groupcounts on a pre-sorted array
 *---------------------------------------------------------------*/
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

/*---------------------------------------------------------------
 * Binary search: find index of 'key' in sorted array, or -1
 *---------------------------------------------------------------*/
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

/*---------------------------------------------------------------
 * Cross-correlation with binary search + OpenMP
 *---------------------------------------------------------------*/
static void cross_correlate_fast(
    const long long* ref_pos, const int* ref_val, int ref_n,
    const long long* tgt_pos, const int* tgt_val, int tgt_n,
    const long long* sft_arr, const int* idx_arr, int num_shifts,
    double* output)
{
#ifdef _OPENMP
    int ncores = omp_get_num_procs();
    int tn = num_shifts < ncores ? num_shifts : ncores;
    if (tn < 1) tn = 1;
#else
    int tn = 1;
#endif

    #pragma omp parallel for schedule(static) num_threads(tn)
    for (int i = 0; i < num_shifts; i++) {
        long long sft = sft_arr[i];
        long long temp = 0;

        for (int j = 0; j < tgt_n; j++) {
            long long key = tgt_pos[j] + sft;
            int idx = binary_find(ref_pos, ref_n, key);
            if (idx >= 0) {
                temp += (long long)ref_val[idx] * tgt_val[j];
            }
        }
        output[idx_arr[i]] = (double)temp;
    }
}

/*---------------------------------------------------------------
 * Compute DFT peak for sub-pixel estimation on ROI_size x ROI_size
 *---------------------------------------------------------------*/
static void dft_peak(const double* ROIcc, int N, double* PX, double* PY) {
    double freq = 2.0 * M_PI / (double)N;

    /* fft_values(0,1) after transpose */
    double re_x = 0.0, im_x = 0.0;
    for (int c = 0; c < N; c++) {
        double s = 0.0;
        for (int j = 0; j < N; j++) {
            s += ROIcc[j * N + c];
        }
        re_x += s * cos(freq * c);
        im_x -= s * sin(freq * c);
    }

    /* fft_values(1,0) after transpose */
    double re_y = 0.0, im_y = 0.0;
    for (int r = 0; r < N; r++) {
        double s = 0.0;
        for (int i = 0; i < N; i++) {
            s += ROIcc[r * N + i];
        }
        re_y += s * cos(freq * r);
        im_y -= s * sin(freq * r);
    }

    double angX = atan2(im_x, re_x);
    double angY = atan2(im_y, re_y);
    if (angX > 0) angX -= 2.0 * M_PI;
    if (angY > 0) angY -= 2.0 * M_PI;

    *PX = (fabs(angX) / (2.0 * M_PI / N) + 1.0) - (N + 1.0) / 2.0;
    *PY = (fabs(angY) / (2.0 * M_PI / N) + 1.0) - (N + 1.0) / 2.0;
}

/*---------------------------------------------------------------
 * MEX gateway
 *---------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 9)
        mexErrMsgIdAndTxt("IntersectionMaxCore:nrhs", "Nine inputs required.");

    double* XList    = mxGetPr(prhs[0]);
    double* YList    = mxGetPr(prhs[1]);
    double* refXList = mxGetPr(prhs[2]);
    double* refYList = mxGetPr(prhs[3]);
    double* fID_d    = mxGetPr(prhs[4]);
    int trackNUM     = (int)mxGetScalar(prhs[5]);
    int trackInterval = (int)mxGetScalar(prhs[6]);
    double imW       = mxGetScalar(prhs[7]);
    double IntersectD = mxGetScalar(prhs[8]);

    int nPoints = (int)mxGetNumberOfElements(prhs[0]);
    int nRef    = (int)mxGetNumberOfElements(prhs[2]);

    double scale = imW / IntersectD;
    int roiR = 3;
    int ROI_size = 2 * roiR + 1;
    int num_shifts = ROI_size * ROI_size;

    /* ============================================================
     * OPTIMIZATION 1: Pre-encode all positions to 1D ONCE
     * ============================================================*/
    std::vector<long long> pList(nPoints);
    std::vector<int> fID(nPoints);
    for (int i = 0; i < nPoints; i++) {
        pList[i] = (long long)round(YList[i] / IntersectD) * (long long)scale
                 + (long long)round(XList[i] / IntersectD);
        fID[i] = (int)fID_d[i];
    }

    /* ============================================================
     * OPTIMIZATION 2: Pre-bin points by segment ONCE
     * Instead of scanning all N points for each segment,
     * sort by frame and use offsets.
     * ============================================================*/

    /* Create index array sorted by frame ID */
    std::vector<int> sortIdx(nPoints);
    for (int i = 0; i < nPoints; i++) sortIdx[i] = i;
    std::sort(sortIdx.begin(), sortIdx.end(),
        [&fID](int a, int b) { return fID[a] < fID[b]; });

    /* Build segment start/end offsets */
    /* seg_start[s] = first index in sortIdx where fID > s*trackInterval
     * seg_end[s]   = last index where fID <= (s+1)*trackInterval */
    std::vector<int> seg_start(trackNUM + 1, nPoints);
    std::vector<int> seg_end(trackNUM + 1, 0);

    {
        int s = 0;
        for (int i = 0; i < nPoints; i++) {
            int f = fID[sortIdx[i]];
            int seg = (f - 1) / trackInterval;  /* 0-based segment */
            if (seg >= trackNUM) seg = trackNUM - 1;

            if (i < seg_start[seg]) seg_start[seg] = i;
            if (i + 1 > seg_end[seg]) seg_end[seg] = i + 1;
        }
    }

    /* ============================================================
     * Encode reference and groupcounts (sorted)
     * ============================================================*/
    std::vector<long long> refEncoded(nRef);
    for (int i = 0; i < nRef; i++) {
        refEncoded[i] = (long long)round(refYList[i] / IntersectD) * (long long)scale
                      + (long long)round(refXList[i] / IntersectD);
    }
    std::sort(refEncoded.begin(), refEncoded.end());

    std::vector<long long> ref_vals;
    std::vector<int> ref_cnts;
    groupcounts_sorted(refEncoded, ref_vals, ref_cnts);
    int ref_n = (int)ref_vals.size();

    /* Build shift and index arrays */
    std::vector<long long> array_sft(num_shifts);
    std::vector<int> array_idx(num_shifts);
    int num = 0;
    for (int r = -roiR; r <= roiR; r++) {
        for (int c = -roiR; c <= roiR; c++) {
            array_sft[num] = (long long)(r * scale + c);
            array_idx[num] = (r + roiR) * ROI_size + (c + roiR);
            num++;
        }
    }

    /* Allocate output */
    plhs[0] = mxCreateDoubleMatrix(1, trackNUM, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, trackNUM, mxREAL);
    double* driftX = mxGetPr(plhs[0]);
    double* driftY = mxGetPr(plhs[1]);

    /* Pre-allocate segment buffers */
    std::vector<long long> seg_encoded;
    std::vector<long long> seg_vals;
    std::vector<int> seg_cnts;
    std::vector<double> img_crr(num_shifts);

    seg_encoded.reserve(nPoints / trackNUM * 2);

    /* ============================================================
     * Main tracking loop
     * ============================================================*/
    double refx = 0.0, refy = 0.0;

    for (int s = 1; s < trackNUM; s++) {
        /* OPTIMIZATION 2 payoff: O(segment_size) instead of O(nPoints) */
        seg_encoded.clear();
        if (seg_start[s] < seg_end[s]) {
            for (int i = seg_start[s]; i < seg_end[s]; i++) {
                int idx = sortIdx[i];
                /* Double-check frame range: fID > s*trackInterval && fID <= (s+1)*trackInterval */
                if (fID[idx] > s * trackInterval && fID[idx] <= (s + 1) * trackInterval) {
                    seg_encoded.push_back(pList[idx]);
                }
            }
        }

        if (seg_encoded.empty()) {
            driftX[s] = (s > 0) ? driftX[s-1] : 0.0;
            driftY[s] = (s > 0) ? driftY[s-1] : 0.0;
            continue;
        }

        /* Sort segment and groupcounts */
        std::sort(seg_encoded.begin(), seg_encoded.end());
        groupcounts_sorted(seg_encoded, seg_vals, seg_cnts);
        int seg_n = (int)seg_vals.size();

        /* Apply coarse shift */
        long long sft = (long long)(round(refy) * scale + round(refx));
        for (int i = 0; i < seg_n; i++) {
            seg_vals[i] += sft;
        }

        /* Cross-correlate with binary search */
        std::fill(img_crr.begin(), img_crr.end(), 0.0);
        cross_correlate_fast(
            ref_vals.data(), ref_cnts.data(), ref_n,
            seg_vals.data(), seg_cnts.data(), seg_n,
            array_sft.data(), array_idx.data(), num_shifts,
            img_crr.data());

        /* DFT peak estimation */
        double PX, PY;
        dft_peak(img_crr.data(), ROI_size, &PX, &PY);

        refx = round(refx) + PX;
        refy = round(refy) + PY;
        driftX[s] = -refx;
        driftY[s] = -refy;
    }
}