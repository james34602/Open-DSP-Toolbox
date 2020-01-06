// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; see the file COPYING. If not, see
// <https://www.gnu.org/licenses/>.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "fht.h"
#include "solveLinearSystem.h"
#include "cpoly.h"
#include "iirToolbox.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
void InitEquationErrorIIR(EquationErrorIIR *iir, int num, int denom, int reqLen)
{
	iir->num = num;
	iir->denom = denom;
	iir->reqLen = reqLen;
	iir->A0xLen = num + denom + 1;
	iir->twoReqLen = reqLen << 1;
	int length = ((reqLen * iir->A0xLen) << 1) + (reqLen * 6) + iir->twoReqLen + iir->twoReqLen * iir->A0xLen + iir->A0xLen;
	iir->memSize = (length + (num + 1) + (denom + 1)) * sizeof(double);
	iir->memoryBuffer = (double*)malloc(iir->memSize);
	memset(iir->memoryBuffer, 0, iir->memSize);
	iir->b = iir->memoryBuffer + length;
	iir->a = iir->b + (num + 1);
}
void EquationErrorIIRFree(EquationErrorIIR *iir)
{
	free(iir->memoryBuffer);
}
void eqnerror(EquationErrorIIR *iir, double *om, double *DReal, double *DImag, double *W, int iter)
{
	if (iter < 1)
		iter = 1;
	int i, j;
	int num = iir->num, denom = iir->denom, reqLen = iir->reqLen, A0xLen = iir->A0xLen, twoReqLen = iir->twoReqLen;
	memset(iir->memoryBuffer, 0, iir->memSize);
	double *A0Re = iir->memoryBuffer;
	double *A0Im = A0Re + (reqLen * A0xLen);
	double *W0 = A0Im + (reqLen * A0xLen);
	memcpy(W0, W, reqLen * sizeof(double));
	double *pvRe = W0 + reqLen;
	double *pvIm = pvRe + reqLen;
	double *sRe = pvIm + reqLen;
	double *sIm = sRe + reqLen;
	for (i = 0; i < reqLen; i++)
	{
		sRe[i] = cos(om[i]);
		sIm[i] = sin(om[i]);
	}
	double *den = sIm + reqLen;
	double *weightedCplxResponse = den + reqLen;
	double *weightedCplxLinearSystem = weightedCplxResponse + twoReqLen;
	double *x = weightedCplxLinearSystem + iir->twoReqLen * iir->A0xLen;
	double *b = iir->b;
	double *a = iir->a;
	a[0] = 1.0;
	for (i = 0; i < reqLen; i++)
	{
		den[i] = 1.0;
		for (j = 0; j < denom; j++)
		{
			double real = cos(om[i] * (double)(j + 1));
			double imag = -sin(om[i] * (double)(j + 1));
			A0Re[i * A0xLen + j] = -(DReal[i] * real - DImag[i] * imag);
			A0Im[i * A0xLen + j] = -(DReal[i] * imag + DImag[i] * real);
		}
		for (j = 0; j < num + 1; j++)
		{
			A0Re[i * A0xLen + denom + j] = cos(om[i] * (double)j);
			A0Im[i * A0xLen + denom + j] = -sin(om[i] * (double)j);
		}
	}
	for (int it = 0; it < iter; it++)
	{
		for (i = 0; i < reqLen; i++)
		{
			W0[i] /= den[i];
			for (j = 0; j < A0xLen; j++)
			{
				weightedCplxLinearSystem[i * A0xLen + j] = A0Re[i * A0xLen + j] * W0[i];
				weightedCplxLinearSystem[(i + reqLen) * A0xLen + j] = A0Im[i * A0xLen + j] * W0[i];
			}
			weightedCplxResponse[i] = DReal[i] * W0[i];
			weightedCplxResponse[i + reqLen] = DImag[i] * W0[i];
		}
		int matSize[2];
		mldivide(weightedCplxLinearSystem, iir->twoReqLen, A0xLen, weightedCplxResponse, iir->twoReqLen, 1, x, matSize);
		for (i = denom; i < A0xLen; i++)
			b[i - denom] = x[i];
		for (i = 1; i < denom + 1; i++)
			a[i] = x[i - 1];
		memset(pvIm, 0, reqLen * sizeof(double));
		for (i = 0; i < reqLen; i++)
			pvRe[i] = 1.0;
		for (i = 0; i < denom; i++)
			for (j = 0; j < reqLen; j++)
			{
				double s_im = sRe[j] * pvIm[j] + sIm[j] * pvRe[j];
				pvRe[j] = sRe[j] * pvRe[j] - sIm[j] * pvIm[j] + a[i + 1];
				pvIm[j] = s_im;
			}
		for (i = 0; i < reqLen; i++)
		{
			double omexpdenomRe = cos(om[i] * (double)denom);
			double imag = sin(om[i] * (double)denom);
			double real = (pvRe[i] * omexpdenomRe + pvIm[i] * imag) / (omexpdenomRe * omexpdenomRe + imag * imag);
			imag = (pvIm[i] * omexpdenomRe - pvRe[i] * imag) / (omexpdenomRe * omexpdenomRe + imag * imag);
			den[i] = sqrt(real * real + imag * imag);
		}
	}
}
// Matlab + Octave style transfer function conversion function
// Stress test proved, shouldn't leak memory/crash.
// All function behaviour is very similar to octave
// zp2sos
void swap(double* a, double* b)
{
	double t = *a;
	*a = *b;
	*b = t;
}
void selectionSort(double arr[], int n)
{
	int i, j, min_idx;
	for (i = 0; i < n - 1; i++)
	{
		min_idx = i;
		for (j = i + 1; j < n; j++)
			if (arr[j] < arr[min_idx])
				min_idx = j;
		swap(&arr[min_idx], &arr[i]);
	}
}
void selectionSortAux(double arr[], double arr2[], int n)
{
	int i, j, min_idx;
	for (i = 0; i < n - 1; i++)
	{
		min_idx = i;
		for (j = i + 1; j < n; j++)
			if (arr[j] < arr[min_idx])
				min_idx = j;
		swap(&arr[min_idx], &arr[i]);
		swap(&arr2[min_idx], &arr2[i]);
	}
}
int cplxpair(double *xRe, double *xIm, unsigned int xLen, double *sortedRe, double *sortedIm)
{
	unsigned int i, j;
	double tol = 100.0 * DBL_EPSILON;
	double *xcRe = (double*)malloc(xLen * sizeof(double));
	double *xcIm = (double*)malloc(xLen * sizeof(double));
	memcpy(xcRe, xRe, xLen * sizeof(double));
	memcpy(xcIm, xIm, xLen * sizeof(double));
	unsigned int *aryIdx = (int*)malloc(xLen * sizeof(int));
	double *tmp1 = (double*)malloc(xLen * sizeof(double));
	unsigned int index = 0;
	while (1) // Odd number of entries remaining
	{
		for (i = 0; i < xLen; i++)
		{
			if (fabs(xcIm[i]) <= tol * sqrt(xcRe[i] * xcRe[i] + xcIm[i] * xcIm[i]))
			{
				aryIdx[index] = i;
				tmp1[index++] = xcRe[i];
			}
		}
		if ((xLen - index) % 2 != 0)
		{
			index = 0;
			tol *= 10.0;
			continue;
		}
		else
			break;
	}
	index = 0;
	for (i = 0; i < xLen; i++)
	{
		if (fabs(xcIm[i]) <= tol * sqrt(xcRe[i] * xcRe[i] + xcIm[i] * xcIm[i]))
		{
			aryIdx[index] = i;
			tmp1[index++] = xcRe[i];
		}
	}
	selectionSort(tmp1, index);
	for (i = xLen - index; i < xLen; i++)
	{
		sortedRe[i] = tmp1[i - xLen + index];
		sortedIm[i] = 0.0;
	}
	unsigned int loop = 0;
	for (i = 0; i < index; i++)
	{
		for (unsigned int idx = 0; idx < xLen + loop; idx++)
		{
			if (idx == aryIdx[i])
			{
				for (j = idx - loop; j < xLen - 1; j++)
				{
					xcRe[j] = xcRe[j + 1];
					xcIm[j] = xcIm[j + 1];
				}
				xLen--;
				loop++;
			}
		}
	}
	if (!xLen)
	{
		free(xcRe);
		free(xcIm);
		free(aryIdx);
		free(tmp1);
		return -2;
	}
	if (xLen % 2 != 0) // Odd number of entries remaining
	{
		printf("Complex numbers can't be paired. ");
		free(xcRe);
		free(xcIm);
		free(aryIdx);
		free(tmp1);
		return -1;
	}
	// Sort complex column-vector xc, based on its real part
	selectionSortAux(xcRe, xcIm, xLen);
	// Check real part pairs to see if imag parts are conjugates
	unsigned int nxt_row = 0; // next row in y for results
	double tmp2[2];
	tol = 100.0 * DBL_EPSILON;
	int previousFail = 0;
	while (xLen)
	{
		unsigned int nn = 0;
		for (i = 0; i < xLen; i++)
			if (fabs(xcRe[i] - xcRe[0]) <= tol * sqrt(xcRe[i] * xcRe[i] + xcIm[i] * xcIm[i]))
				aryIdx[nn++] = i;
		if (nn <= 1 || nn > 2)
		{
			if (tol > 1e-5)
				break; // Simply no complex numbers pair
			tol *= 10.0;
			printf("Complex numbers can't be paired, continue with larger tolerance\n");
			previousFail = 1;
			continue;
		}
		else
		{
			tol = 100.0 * DBL_EPSILON;
			previousFail = 0;
		}
		for (i = 0; i < nn; i++)
		{
			tmp1[i] = xcIm[i];
			tmp2[i] = xcRe[i];
		}
		selectionSortAux(tmp1, tmp2, nn);
		sortedRe[nxt_row] = tmp2[1];
		sortedIm[nxt_row] = -tmp1[1];
		sortedRe[nxt_row + 1] = tmp2[1];
		sortedIm[nxt_row + 1] = tmp1[1];
		nxt_row += nn;
		loop = 0;
		for (i = 0; i < nn; i++)
		{
			for (unsigned int idx = 0; idx < xLen + loop; idx++)
			{
				if (idx == aryIdx[i])
				{
					for (j = idx - loop; j < xLen - 1; j++)
					{
						xcRe[j] = xcRe[j + 1];
						xcIm[j] = xcIm[j + 1];
					}
					xLen--;
					loop++;
				}
			}
		}
	}
	free(xcRe);
	free(xcIm);
	free(aryIdx);
	free(tmp1);
	return 0;
}
int zp2sos(double *zRe, double *zIm, unsigned int zLen, double *pRe, double *pIm, unsigned int pLen, double *sos)
{
	unsigned int i;
	const double thresh = 100 * DBL_EPSILON;
	unsigned int nzc = 0, nzr = 0, npc = 0, npr = 0;
	double *zcpRe = (double*)malloc(zLen * sizeof(double));
	double *zcpIm = (double*)malloc(zLen * sizeof(double));
	unsigned int nzrsec = 0, idx;
	if (zLen)
	{
		cplxpair(zRe, zIm, zLen, zcpRe, zcpIm); // sort complex pairs, real roots at end
		idx = zLen - 1;
		while ((idx + 1) && fabs(zcpIm[idx]) < thresh) // determine no.of real values
		{
			nzrsec = nzrsec + 1;
			idx = idx - 1;
		}
	}
	unsigned int nzsect2 = zLen - nzrsec;
	if (nzsect2 % 2 != 0)
	{
		printf("Odd number of zeros!");
		free(zcpRe);
		free(zcpIm);
		return -1;
	}
	nzc = nzsect2 >> 1;
	double *zcRe = (double*)malloc(nzc * sizeof(double));
	double *zcIm = (double*)malloc(nzc * sizeof(double));
	idx = 0;
	for (i = 0; i < nzsect2; i++)
	{
		if ((i + 1) % 2 == 0)
		{
			zcRe[idx] = zcpRe[i];
			zcIm[idx++] = zcpIm[i];
		}
	}
	nzr = zLen - nzsect2;
	double *zr = (double*)malloc((nzr + 1) * sizeof(double));
	for (i = nzsect2; i < zLen; i++)
		zr[i - nzsect2] = zcpRe[i];
	free(zcpRe);
	free(zcpIm);
	zcpRe = (double*)malloc(pLen * sizeof(double));
	zcpIm = (double*)malloc(pLen * sizeof(double));
	nzrsec = 0;
	if (pLen)
	{
		cplxpair(pRe, pIm, pLen, zcpRe, zcpIm); // sort complex pairs, real roots at end
		idx = pLen - 1;
		while ((idx + 1) && fabs(zcpIm[idx]) < thresh) // determine no.of real values
		{
			nzrsec = nzrsec + 1;
			idx = idx - 1;
		}
	}
	nzsect2 = pLen - nzrsec;
	if (nzsect2 % 2 != 0)
	{
		printf("Odd number of zeros!");
		free(zcpRe);
		free(zcpIm);
		free(zcRe);
		free(zcIm);
		free(zr);
		return -1;
	}
	npc = nzsect2 >> 1;
	double *pcRe = (double*)malloc(npc * sizeof(double));
	double *pcIm = (double*)malloc(npc * sizeof(double));
	idx = 0;
	for (i = 0; i < nzsect2; i++)
	{
		if ((i + 1) % 2 == 0)
		{
			pcRe[idx] = zcpRe[i];
			pcIm[idx++] = zcpIm[i];
		}
	}
	npr = pLen - nzsect2;
	double *pr = (double*)malloc((npr + 1) * sizeof(double));
	for (i = nzsect2; i < pLen; i++)
		pr[i - nzsect2] = zcpRe[i];
	free(zcpRe);
	free(zcpIm);

	// Pair up real zeros:
	double *zrms = 0, *zrp = 0;
	if (nzr)
	{
		if (nzr % 2 != 0)
		{
			nzr++;
			zr[nzr - 1] = 0.0;
		}
		nzrsec = nzr >> 1;
		zrms = (double*)malloc(nzrsec * sizeof(double));
		zrp = (double*)malloc(nzrsec * sizeof(double));
		idx = 0;
		for (i = 0; i < nzr; i++)
		{
			if ((i + 1) % 2 != 0)
			{
				zrms[idx] = -zr[i] - zr[i + 1];
				zrp[idx++] = zr[i] * zr[i + 1];
			}
		}
	}
	else
		nzrsec = 0;

	// Pair up real poles:
	unsigned int nprsec;
	double *prms = 0, *prp = 0;
	if (npr)
	{
		if (npr % 2 != 0)
		{
			npr++;
			pr[npr - 1] = 0.0;
		}
		nprsec = npr >> 1;
		prms = (double*)malloc(nprsec * sizeof(double));
		prp = (double*)malloc(nprsec * sizeof(double));
		idx = 0;
		for (i = 0; i < npr; i++)
		{
			if ((i + 1) % 2 != 0)
			{
				prms[idx] = -pr[i] - pr[i + 1];
				prp[idx++] = pr[i] * pr[i + 1];
			}
		}
	}
	else
		nprsec = 0;
	unsigned int nzrl = nzc + nzrsec; // index of last real zero section
	unsigned int nprl = npc + nprsec; // index of last real pole section
	unsigned int nsecs = max(nzrl, nprl);
	// Convert complex zeros and poles to real 2nd-order section form:
	for (i = 0; i < nsecs; i++)
	{
		sos[i * 6] = sos[i * 6 + 3] = 1.0;
		if (i < nzc) // lay down a complex zero pair:
		{
			sos[i * 6 + 1] = -2.0 * zcRe[i];
			sos[i * 6 + 2] = zcRe[i] * zcRe[i] + zcIm[i] * zcIm[i];
		}
		else if (i < nzrl) // lay down a pair of real zeros:
		{
			sos[i * 6 + 1] = zrms[i - nzc];
			sos[i * 6 + 2] = zrp[i - nzc];
		}
		if (i < npc) // lay down a complex pole pair:
		{
			sos[i * 6 + 4] = -2.0 * pcRe[i];
			sos[i * 6 + 5] = pcRe[i] * pcRe[i] + pcIm[i] * pcIm[i];
		}
		else if (i < nprl) // lay down a pair of real poles:
		{
			sos[i * 6 + 4] = prms[i - npc];
			sos[i * 6 + 5] = prp[i - npc];
		}
	}
	free(zcRe);
	free(zcIm);
	free(zr);
	free(pcRe);
	free(pcIm);
	free(pr);
	if (zrms)
		free(zrms);
	if (zrp)
		free(zrp);
	if (prms)
		free(prms);
	if (prp)
		free(prp);
	return nsecs;
}
int tf2sos(double *b, unsigned int bLen, double *a, unsigned int aLen, double **sos)
{
	// % Find Poles and Zeros
	if (!aLen)
	{
		printf("Denominator cannot be empty\n");
		return -1;
	}
	if (a[0] == 0.0)
	{
		printf("Denominator cannot be zero\n");
		return -1;
	}
	unsigned int i;
	double *bpolyRe = (double*)malloc(bLen * sizeof(double));
	memcpy(bpolyRe, b, bLen * sizeof(double));
	double *bpolyIm = (double*)malloc(bLen * sizeof(double));
	memset(bpolyIm, 0, bLen * sizeof(double));
	double *zRe = (double*)malloc(bLen * sizeof(double));
	double *zIm = (double*)malloc(bLen * sizeof(double));
	int zeroNumRoots = cpoly(bpolyRe, bpolyIm, bLen - 1, zRe, zIm);
	free(bpolyRe);
	free(bpolyIm);
	if (zeroNumRoots < 0)
		zeroNumRoots = 0;
	double *apolyRe = (double*)malloc(aLen * sizeof(double));
	memcpy(apolyRe, a, aLen * sizeof(double));
	double *apolyIm = (double*)malloc(aLen * sizeof(double));
	memset(apolyIm, 0, aLen * sizeof(double));
	double *pRe = (double*)malloc(aLen * sizeof(double));
	double *pIm = (double*)malloc(aLen * sizeof(double));
	int poleNumRoots = cpoly(apolyRe, apolyIm, aLen - 1, pRe, pIm);
	free(apolyRe);
	free(apolyIm);
	if (poleNumRoots < 0)
		poleNumRoots = 0;
	unsigned int myNSec = (unsigned int)ceil(max(zeroNumRoots, poleNumRoots) * 0.5);
	*(sos) = (double*)malloc(myNSec * 6 * sizeof(double));
	memset(*(sos), 0, myNSec * 6 * sizeof(double));
	double firstNonZero = 0.0;
	for (i = 0; i < bLen; i++)
	{
		if (b[i] != 0.0)
		{
			firstNonZero = b[i];
			break;
		}
	}
	double k = firstNonZero / a[0];
	int numSections = zp2sos(zRe, zIm, zeroNumRoots, pRe, pIm, poleNumRoots, *(sos));
	free(zRe);
	free(zIm);
	free(pRe);
	free(pIm);
	(*sos)[0] *= k;
	(*sos)[1] *= k;
	(*sos)[2] *= k;
	return numSections;
}
inline void diff(double *y, double *f, int sz)
{
	--sz;
	for (int i = 0; i < sz; i++)
		f[i] = y[i + 1] - y[i];
}
inline int isneg(double *y, int sz)
{
	for (int i = 0; i < sz; i++)
		if (y[i] < 0) return 1;
	return 0;
}
// Design arbitrary response IIR filter
// Function interpolate user provided frequency and gain vector into larger grid
// function perform minimum phase transform to interpolated grid
// arbitrary response filter designing algorithm will be applied
// int gridLen: Interpolation grid points
// double *ff: Frequency vector [0 ... 1]
// double *aa: Gain vector in dB
// int gVLen: Frequency or Gain vector length, length(ff) must == length(aa), unexpected error otherwise
// double *b: Returning array that contain numerator
// int M: Numerator orders, recommended range [2 - 23]
// double *a: Returning array that contain denominator 
// int N: Denominator orders, recommended range [2 - 23]
// int iterEqnErr: Minimization iteration counts, recommended range[1 - 2]
// Returning error code: All error code smaller than 1 mean error detected, no valuable filter have been generated
int designMinimumPhaseArbIIR(int gridLen, double *ff, double *aa, int gVLen, double *b, int M, double *a, int N, int iterEqnErr)
{
	int npt, lap, npt2;
	// Convert gain to linear scale
	for (lap = 0; lap < gVLen; lap++)
		aa[lap] = pow(10.0, aa[lap] / 20.0);
	if (gridLen < 1024)
		npt = 512;
	else
		npt = (int)pow(2.0, ceil(log((double)gridLen) / log(2.0)));
	lap = (int)(npt / 25);
	if (ff[0] != 0.0 || ff[gVLen - 1] != 1.0)
	{
		printf("The first frequency must be 0 and the last 1");
		return 0;
	}
	// Interpolate breakpoints onto large grid
	int nint = gVLen - 1;
	double *df = (double*)malloc(sizeof(double)*(gVLen - 1));
	diff(ff, df, gVLen);
	if (isneg(df, gVLen - 1))
	{
		printf("Frequencies must be non-decreasing");
		return -1;
	}
	npt = npt; // Length of [dc 1 2 ... nyquist] frequencies.
	npt2 = npt * 2;
	int nb = 0;
	int ne = 0;
	double inc;
	double *H = (double*)malloc(sizeof(double)*npt);
	H[0] = aa[0];
	int i;
	for (i = 0; i < nint; i++)
	{
		if (df[i] == 0)
		{
			nb = nb - lap / 2;
			ne = nb + lap;
		}
		else
			ne = (int)(ff[i + 1] * npt) - 1;
		if (nb < 0 || ne > npt)
		{
			printf("Too abrupt an amplitude change near end of frequency interval");
			return -2;
		}
		int j;
		for (j = nb; j <= ne; j++)
		{
			if (nb == ne)
				inc = 0;
			else
				inc = (double)(j - nb) / (double)(ne - nb);
			H[j] = inc * aa[i + 1] + (1.0 - inc)*aa[i];
		}
		nb = ne + 1;
	}
	free(df);
	const double threshdB = -140.0;
	const double threshold = pow(10.0, threshdB / 20.0);
	double logThreshold = log(threshold);
	unsigned int *mBitRev = (unsigned int*)malloc(npt2 * sizeof(unsigned int));
	double *mSineTab = (double*)malloc(npt2 * sizeof(double));
	double *timeData = (double*)malloc(npt2 * sizeof(double));
	double *freqData = (double*)malloc(npt2 * sizeof(double));
	LLbitReversalTbl(mBitRev, npt2);
	LLsinHalfTbl(mSineTab, npt2);
	for (i = 0; i < npt; i++)
	{
		double gain = H[i];
		// Part of minimum phase spectrum --- cepstrum calculation
		if (gain < threshold)
			gain = logThreshold;
		else
			gain = log(gain);
		freqData[mBitRev[i]] = gain;
		freqData[mBitRev[npt2 - 1 - i]] = gain;
	}
	free(H);
	LLdiscreteHartley(freqData, npt2, mSineTab);
	double gain = 1.0 / ((double)npt2);
	timeData[0] = freqData[0] * gain;
	timeData[mBitRev[npt]] = freqData[npt] * gain;
	for (i = 1; i < npt; i++)
	{
		timeData[mBitRev[i]] = (freqData[i] + freqData[npt2 - i]) * gain;
		timeData[mBitRev[npt2 - i]] = 0.0f;
	}
	LLdiscreteHartley(timeData, npt2, mSineTab);
	int gridLength = npt + 1;
	double *om = (double*)malloc(gridLength * sizeof(double));
	double *DRe = (double*)malloc(gridLength * sizeof(double));
	double *DIm = (double*)malloc(gridLength * sizeof(double));
	double *W = (double*)malloc(gridLength * sizeof(double));
	DRe[0] = exp(timeData[0]);
	DIm[0] = 0.0;
	om[0] = 0.0;
	W[0] = 1.0;
	for (i = 1; i < npt + 1; i++)
	{
		om[i] = ((double)i * (1.0 / (double)npt)) * M_PI;
		W[i] = 1.0;
		double re = (timeData[i] + timeData[npt2 - i]) * 0.5;
		double im = (timeData[i] - timeData[npt2 - i]) * 0.5;
		double eR = exp(re);
		DRe[i] = eR * cos(im);
		DIm[i] = -eR * sin(im);
	}
	EquationErrorIIR initSolution;
	InitEquationErrorIIR(&initSolution, M, N, gridLength);
	eqnerror(&initSolution, om, DRe, DIm, W, iterEqnErr);
	memcpy(b, initSolution.b, (M + 1) * sizeof(double));
	memcpy(a, initSolution.a, (N + 1) * sizeof(double));
	EquationErrorIIRFree(&initSolution);
	free(mBitRev);
	free(mSineTab);
	free(timeData);
	free(freqData);
	free(om);
	free(DRe);
	free(DIm);
	free(W);
	return 1;
}