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
#include <time.h>
#include <vld.h>
#include "cpoly.h"
void xzlarf(int cols1, int rows1, int iv0, double tau, double C_data[], int ic0, int ldc, double work_data[])
{
	int lastv;
	int lastc;
	int i;
	char exitg2;
	int jy;
	int i3;
	int ia;
	int ix;
	int exitg1;
	double c;
	int i4;
	int ijA;
	if (tau != 0.0)
	{
		lastv = cols1;
		i = iv0 + cols1;
		while ((lastv > 0) && (C_data[i - 2] == 0.0))
		{
			lastv--;
			i--;
		}
		lastc = rows1 - 1;
		exitg2 = 0;
		while ((!exitg2) && (lastc + 1 > 0))
		{
			i = ic0 + lastc * ldc;
			ia = i;
			do
			{
				exitg1 = 0;
				if (ia <= (i + lastv) - 1)
				{
					if (C_data[ia - 1] != 0.0)
						exitg1 = 1;
					else
						ia++;
				}
				else
				{
					lastc--;
					exitg1 = 2;
				}
			} while (exitg1 == 0);
			if (exitg1 == 1)
				exitg2 = 1;
		}
	}
	else
	{
		lastv = 0;
		lastc = -1;
	}
	if (lastv > 0)
	{
		if (lastc + 1 != 0)
		{
			if (0 <= lastc)
				memset(&work_data[0], 0, (lastc + 1) * sizeof(double));
			i = 0;
			i3 = ic0 + ldc * lastc;
			for (jy = ic0; ldc < 0 ? jy >= i3 : jy <= i3; jy += ldc)
			{
				ix = iv0;
				c = 0.0;
				i4 = (jy + lastv) - 1;
				for (ia = jy; ia <= i4; ia++)
				{
					c += C_data[ia - 1] * C_data[ix - 1];
					ix++;
				}
				work_data[i] += c;
				i++;
			}
		}
		if (-tau != 0.0)
		{
			i = ic0 - 1;
			jy = 0;
			for (ia = 0; ia <= lastc; ia++)
			{
				if (work_data[jy] != 0.0)
				{
					c = work_data[jy] * -tau;
					ix = iv0;
					i3 = i + 1;
					i4 = lastv + i;
					for (ijA = i3; ijA <= i4; ijA++)
					{
						C_data[ijA - 1] += C_data[ix - 1] * c;
						ix++;
					}
				}
				jy++;
				i += ldc;
			}
		}
	}
}
void xswap(int rows1, double x_data[], int ix0, int iy0)
{
	int ix;
	int iy;
	int k;
	double temp;
	ix = ix0 - 1;
	iy = iy0 - 1;
	for (k = 0; k < rows1; k++)
	{
		temp = x_data[ix];
		x_data[ix] = x_data[iy];
		x_data[iy] = temp;
		ix++;
		iy++;
	}
}
double xnrm2(int rows1, const double x_data[], int ix0)
{
	double y;
	double scale;
	int kend;
	int k;
	double absxk;
	double t;
	y = 0.0;
	if (rows1 >= 1)
	{
		if (rows1 == 1)
			y = fabs(x_data[ix0 - 1]);
		else
		{
			scale = DBL_MIN;
			kend = (ix0 + rows1) - 1;
			for (k = ix0; k <= kend; k++)
			{
				absxk = fabs(x_data[k - 1]);
				if (absxk > scale)
				{
					t = scale / absxk;
					y = 1.0 + y * t * t;
					scale = absxk;
				}
				else
				{
					t = absxk / scale;
					y += t * t;
				}
			}
			y = scale * sqrt(y);
		}
	}
	return y;
}
static double rt_hypotd(double u0, double u1)
{
	double y;
	double a;
	double b;
	a = fabs(u0);
	b = fabs(u1);
	if (a < b)
	{
		a /= b;
		y = b * sqrt(a * a + 1.0);
	}
	else if (a > b)
	{
		b /= a;
		y = a * sqrt(b * b + 1.0);
	}
	else
		y = a * sqrt(2.0);
	return y;
}
void xgeqp3(double A[], const int rows1, const int cols1, double tau_data[], int jpvt_data[], double *work_data)
{
	int b_n;
	int mn;
	int yk;
	int k;
	double *vn1_data = work_data + rows1;
	double *vn2_data = vn1_data + rows1;
	int i;
	double smax;
	int i_i;
	int nmi;
	int mmi;
	double s;
	double beta1;
	int knt;
	b_n = cols1;
	mn = rows1;
	if (b_n < mn)
		mn = b_n;
	jpvt_data[0] = 1;
	yk = 1;
	for (k = 2; k <= rows1; k++)
		jpvt_data[k - 1] = ++yk;
	if (0 <= rows1 - 1)
		memset(&work_data[0], 0, rows1 * (int)sizeof(double));
	k = 1;
	for (yk = 0; yk < rows1; yk++)
	{
		smax = xnrm2(cols1, A, k);
		vn1_data[yk] = smax;
		vn2_data[yk] = smax;
		k += cols1;
	}
	for (i = 0; i < mn; i++)
	{
		i_i = i + i * cols1;
		nmi = rows1 - i;
		mmi = (cols1 - i) - 1;
		if (nmi < 1)
			b_n = 0;
		else
		{
			b_n = 1;
			if (nmi > 1)
			{
				yk = i;
				smax = fabs(vn1_data[i]);
				for (k = 2; k <= nmi; k++)
				{
					yk++;
					s = fabs(vn1_data[yk]);
					if (s > smax)
					{
						b_n = k;
						smax = s;
					}
				}
			}
		}
		b_n = (i + b_n) - 1;
		if (b_n + 1 != i + 1)
		{
			xswap(cols1, A, 1 + cols1 * b_n, 1 + cols1 * i);
			yk = jpvt_data[b_n];
			jpvt_data[b_n] = jpvt_data[i];
			jpvt_data[i] = yk;
			vn1_data[b_n] = vn1_data[i];
			vn2_data[b_n] = vn2_data[i];
		}
		if (i + 1 < cols1)
		{
			s = A[i_i];
			b_n = i_i + 2;
			tau_data[i] = 0.0;
			if (1 + mmi > 0)
			{
				smax = xnrm2(mmi, A, i_i + 2);
				if (smax != 0.0)
				{
					beta1 = rt_hypotd(A[i_i], smax);
					if (A[i_i] >= 0.0)
						beta1 = -beta1;
					if (fabs(beta1) < DBL_MIN)
					{
						knt = -1;
						do
						{
							knt++;
							for (k = b_n; k <= i_i + mmi + 1; k++)
								A[k - 1] *= DBL_MAX;
							beta1 *= DBL_MAX;
							s *= DBL_MAX;
						} while (!(fabs(beta1) >= DBL_MIN));
						beta1 = rt_hypotd(s, xnrm2(mmi, A, i_i + 2));
						if (s >= 0.0)
							beta1 = -beta1;
						tau_data[i] = (beta1 - s) / beta1;
						smax = 1.0 / (s - beta1);
						for (k = b_n; k <= i_i + mmi + 1; k++)
							A[k - 1] *= smax;
						for (k = 0; k <= knt; k++)
							beta1 *= DBL_MIN;
						s = beta1;
					}
					else
					{
						tau_data[i] = (beta1 - A[i_i]) / beta1;
						smax = 1.0 / (A[i_i] - beta1);
						for (k = b_n; k <= i_i + mmi + 1; k++)
							A[k - 1] *= smax;
						s = beta1;
					}
				}
			}
			A[i_i] = s;
		}
		else
			tau_data[i] = 0.0;
		if (i + 1 < rows1)
		{
			s = A[i_i];
			A[i_i] = 1.0;
			xzlarf(1 + mmi, nmi - 1, i_i + 1, tau_data[i], A, (i + (i + 1) * cols1) + 1, cols1, work_data);
			A[i_i] = s;
		}
		for (yk = i + 2; yk <= rows1; yk++)
		{
			smax = vn1_data[yk - 1];
			if (smax != 0.0)
			{
				b_n = i + cols1 * (yk - 1);
				s = fabs(A[b_n]) / smax;
				s = 1.0 - s * s;
				if (s < 0.0)
					s = 0.0;
				beta1 = smax / vn2_data[yk - 1];
				beta1 = s * (beta1 * beta1);
				if (beta1 <= DBL_EPSILON)
				{
					if (i + 1 < cols1)
					{
						smax = xnrm2(mmi, A, b_n + 2);
						vn1_data[yk - 1] = smax;
						vn2_data[yk - 1] = smax;
					}
					else
					{
						vn1_data[yk - 1] = 0.0;
						vn2_data[yk - 1] = 0.0;
					}
				}
				else
					vn1_data[yk - 1] = smax * sqrt(s);
			}
		}
	}
}
void xgetrf(const int m, double A[], int ipiv_data[])
{
	int yk, jA, u0, j, mmj, b_tmp, jp1j, ix, i2;
	double smax, s;
	ipiv_data[0] = 1;
	yk = 1;
	for (jA = 2; jA <= m; jA++)
		ipiv_data[jA - 1] = ++yk;
	for (j = 0; j < m - 1; j++)
	{
		mmj = m - j;
		b_tmp = j * (m + 1);
		jp1j = b_tmp + 2;
		yk = 0;
		ix = b_tmp;
		smax = fabs(A[b_tmp]);
		for (jA = 2; jA <= mmj; jA++)
		{
			s = fabs(A[++ix]);
			if (s > smax)
			{
				yk = jA - 1;
				smax = s;
			}
		}
		if (A[b_tmp + yk] != 0.0)
		{
			if (yk != 0)
			{
				yk += j;
				ipiv_data[j] = yk + 1;
				ix = j;
				for (jA = 0; jA < m; jA++)
				{
					smax = A[ix];
					A[ix] = A[yk];
					A[yk] = smax;
					ix += m;
					yk += m;
				}
			}
			for (yk = jp1j; yk <= b_tmp + mmj; yk++)
				A[yk - 1] /= A[b_tmp];
		}
		yk = b_tmp + m;
		jA = (b_tmp + m) + 1;
		for (jp1j = 0; jp1j <= m - j - 2; jp1j++)
		{
			smax = A[yk];
			if (A[yk] != 0.0)
			{
				ix = b_tmp + 1;
				i2 = mmj + jA;
				for (u0 = jA + 1; u0 < i2; u0++)
					A[u0 - 1] += A[ix++] * -smax;
			}
			yk += m;
			jA += m;
		}
	}
}
void mldivide(const double A[], const int rows1, const int cols1, const double b[], const int rows2, const int cols2, double Y[], int Y_size[2])
{
	if (rows1 != rows2)
	{
		printf("Matrix dimensions must agree.\n");
		return;
	}
	int minmn, i0, maxmn, m, b_nb, j, mn, k, i;
	double tol;
	size_t bufferSize = ((rows2 * cols2) + (rows1 * cols1) + (cols1 * cols2) + (cols1 * 3) + min(rows1, cols1)) * sizeof(double) + max(rows1, cols1) * sizeof(int);
	char *workingBuffer = (char*)malloc(bufferSize * sizeof(char));
	double *B_tmp_data = (double*)workingBuffer;
	double *b_A_data = B_tmp_data + rows2 * cols2;
	double *b_Y_data = b_A_data + rows1 * cols1;
	double *qrWorkBuf = b_Y_data + cols1 * cols2;
	double *tau_data = qrWorkBuf + cols1 * 3;
	int *jpvt_data = (int*)(tau_data + min(rows1, cols1));
	if (rows1 == cols1)
	{
		for (i = 0; i < cols1; i++)
			for (j = 0; j < rows1; j++)
				b_A_data[j + rows1 * i] = A[i + cols1 * j];
		for (i = 0; i < cols2; i++)
			for (j = 0; j < rows2; j++)
				b_Y_data[j + rows2 * i] = b[i + cols2 * j];
		xgetrf(cols1, b_A_data, jpvt_data);
		for (maxmn = 0; maxmn <= cols1 - 2; maxmn++)
		{
			if (jpvt_data[maxmn] != maxmn + 1)
			{
				for (m = 0; m <= cols2 - 1; m++)
				{
					minmn = rows2 * m;
					b_nb = maxmn + minmn;
					tol = b_Y_data[b_nb];
					mn = jpvt_data[maxmn] + minmn - 1;
					b_Y_data[b_nb] = b_Y_data[mn];
					b_Y_data[mn] = tol;
				}
			}
		}
		for (j = 0; j <= cols2 - 1; j++)
		{
			minmn = cols1 * j - 1;
			for (k = 0; k < cols1; k++)
			{
				if (b_Y_data[k + minmn + 1] != 0.0)
					for (i = k + 2; i <= cols1; i++)
						b_Y_data[i + minmn] -= b_Y_data[k + minmn + 1] * b_A_data[(i + cols1 * k) - 1];
			}
		}
		for (j = 0; j <= cols2 - 1; j++)
		{
			minmn = cols1 * j - 1;
			for (k = cols1; k >= 1; k--)
			{
				maxmn = cols1 * (k - 1) - 1;
				i0 = k + minmn;
				if (b_Y_data[i0] != 0.0)
				{
					b_Y_data[i0] /= b_A_data[k + maxmn];
					for (i = 0; i <= k - 2; i++)
						b_Y_data[i + minmn + 1] -= b_Y_data[i0] * b_A_data[i + maxmn + 1];
				}
			}
		}
		Y_size[1] = cols2;
		Y_size[0] = rows2;
		for (i = 0; i < rows2; i++)
			for (j = 0; j < cols2; j++)
				Y[j + cols2 * i] = b_Y_data[i + rows2 * j];
	}
	else
	{
		for (i = 0; i < cols1; i++)
			for (j = 0; j < rows1; j++)
				b_A_data[j + rows1 * i] = A[i + cols1 * j];
		for (i = 0; i < cols2; i++)
			for (j = 0; j < rows2; j++)
				B_tmp_data[j + rows2 * i] = b[i + cols2 * j];
		xgeqp3(b_A_data, cols1, rows1, tau_data, jpvt_data, qrWorkBuf);
		if (rows1 < cols1)
		{
			minmn = rows1;
			maxmn = cols1;
		}
		else
		{
			minmn = cols1;
			maxmn = rows1;
		}
		int rankR = 0;
		if (minmn > 0)
		{
			tol = DBL_EPSILON * (double)maxmn * fabs(b_A_data[0]);
			while ((rankR < minmn) && (fabs(b_A_data[rankR + rows1 * rankR]) > tol))
				rankR++;
		}
		memset(b_Y_data, 0, cols1 * cols2 * sizeof(double));
		m = rows1;
		b_nb = cols2;
		minmn = rows1;
		mn = cols1;
		if (minmn < mn)
			mn = minmn;
		for (j = 0; j < mn; j++)
		{
			if (tau_data[j] != 0.0)
			{
				for (k = 0; k < b_nb; k++)
				{
					maxmn = rows2 * k;
					minmn = j + maxmn;
					tol = B_tmp_data[minmn];
					i0 = j + 2;
					for (i = i0; i <= m; i++)
						tol += b_A_data[(i + rows1 * j) - 1] * B_tmp_data[(i + maxmn) - 1];
					tol *= tau_data[j];
					if (tol != 0.0)
					{
						B_tmp_data[minmn] -= tol;
						for (i = i0; i <= m; i++)
							B_tmp_data[(i + maxmn) - 1] -= b_A_data[(i + rows1 * j) - 1] * tol;
					}
				}
			}
		}
		for (k = 0; k < cols2; k++)
		{
			for (i = 0; i < rankR; i++)
				b_Y_data[(jpvt_data[i] + cols1 * k) - 1] = B_tmp_data[i + rows2 * k];
			for (j = rankR; j >= 1; j--)
			{
				mn = (jpvt_data[j - 1] + cols1 * k) - 1;
				minmn = rows1 * (j - 1);
				b_Y_data[mn] /= b_A_data[(j + minmn) - 1];
				for (i = 0; i <= j - 2; i++)
				{
					maxmn = (jpvt_data[i] + cols1 * k) - 1;
					b_Y_data[maxmn] -= b_Y_data[mn] * b_A_data[i + minmn];
				}
			}
		}
		Y_size[1] = cols2;
		Y_size[0] = cols1;
		for (i = 0; i < cols1; i++)
			for (j = 0; j < cols2; j++)
				Y[j + cols2 * i] = b_Y_data[i + cols1 * j];
	}
	free(workingBuffer);
}
void linspace(double a, double b, const int n, double *y)
{
	int i;
	double d = (b - a) / (double)(n - 1);
	for (i = 0; i < n; i++)
		y[i] = a + i * d;
}
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef struct
{
	int num, denom, reqLen, A0xLen, twoReqLen;
	size_t memSize;
	double *memoryBuffer;
	double *b, *a;
} EquationErrorIIR;
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
double b1[26] = { 1.0012013530872574, -0.069720018620065285, -0.64673918227041816, -0.55757536000231844, 0.64295878194519718, 0.0094022197292261238, -0.18633856599165166, -0.4691258318656748, 0.4442598323878294, 0.072686761493298849, 0.13409918683169322, -0.12474411858102159, -0.056469100873046407, -0.22006840993047824, 0.26453367463948452, 0.26681588201565581, -0.067606603898450218, -0.32955314592250551, -0.044882442980331214, 0.14387132477091111, 0.019833889782182604, 0.028278334405851192, -0.045872655616353675, 0.0030224372295484503, -0.010703554301967084, -0.01981380052058232 };
double a1[21] = { 1, -0.090427174787398881, -0.47635858646483303, -0.79639935556516706, 0.62662985386978876, -0.019252273869171402, 0.373352197323154, -0.58121888934590016, 0.27218273634698931, -0.37198496466002778, 0.39903095219475698, -0.17678573181944943, 0.3548203543708196, -0.37239762221395434, 0.18161476755467335, -0.073080294381017097, 0.18146689079210729, -0.27364314013427843, 0.040104152008789251, -0.012082609399677246, 0.067263997327562164 };
double b2[26] = { 0.038517370142648298, 0.0010386859687012873, -0.086098050068718457, 0.085348772667813211, -0.042117970125450793, -0.020679539631744896, 0.039125127393445806, 0.014664702676612512, -0.018297543100421525, 0.00016705554654026047, 0.0087178531905720818, -0.018730936858442637, 0.046656201120630895, -0.072544732817444543, 0.024386405662560805, 0.049948360456525134, -0.063828258123357837, 0.041860387613236426, -0.0018425830278492841, 0.01337618002356265, -0.060683481429640203, 0.062629730807545816, -0.035835930065334823, 0.022974602720135216, -0.014730817467697822, 0.0063606915280241668 };
double a2[21] = { 1, -2.6440838619465326, 3.0666426788538801, -2.7142888887443455, 2.4389834313374412, -2.3154425986172842, 2.7314555000350795, -2.7124702506075016, 1.8813573943472346, -1.1376376277240625, 0.56964901530213952, -0.0782985386319653, -0.25201252932238738, 0.31903755140407497, -0.17444732353420295, 0.19097011797437338, -0.31444819363505749, 0.3308160548842019, -0.1889454622592838, 0.0030697514846655016, 0.034858537587783761 };
double b3[12] = { 0.696176763775895,0.555547146510916,0.04300074335326,0.16659501834806,0.586677501045686,0.509389108644681,0.999366940130689,0.978550763831388,0.764782787743648,0.234238613392499,0.736905912661616,0.666684634980556 };
double a3[19] = { 1,0.248893139221197,0.342181213844591,0.219319017669164,0.362032698557814,0.00710355493812713,0.234346093282161,0.112449449721568,0.297554824174743,0.0863682294423248,0.0203849686160456,0.112731654169387,0.053411714901486,0.182745972557301,0.0891860729537674,0.221308153546729,0.0284603569042661,0.0389974631641958,0.129727551826713 };
double b4[4] = { -0.3, 0.2, 0.01, -2.15651 };
double a4[4] = { -7.0, 0.0, 0.0, 0.0 };
// tf2sos tests
/*int main()
{
	int i, j;
	double *sos;
	int numSecs = tf2sos(b1, 26, a1, 21, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	free(sos);
	printf("\n");
	numSecs = tf2sos(b2, 26, a2, 21, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	free(sos);
	printf("\n");
	numSecs = tf2sos(b3, 12, a3, 19, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	free(sos);
	printf("\n");
	numSecs = tf2sos(b4, 4, a4, 4, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	free(sos);
	return 0;
}*/
// Wideband passband linear phase low pass filter
int main()
{
	int i;
	const double tau = 18.0;
	const int M = 14;
	const int N = 11;
	const int gridLength = 88;
	int passband = 40;
	int leftOver = gridLength - passband;
	double *om = (double*)malloc(gridLength * sizeof(double));
	double *DRe = (double*)malloc(gridLength * sizeof(double));
	double *DIm = (double*)malloc(gridLength * sizeof(double));
	memset(DRe, 0, gridLength * sizeof(double));
	memset(DIm, 0, gridLength * sizeof(double));
	double *W = (double*)malloc(gridLength * sizeof(double));
	linspace(0.0, 0.4, passband, om);
	linspace(0.56, 1.0, leftOver, om + passband);
	for (i = 0; i < gridLength; i++)
		om[i] *= M_PI;
	for (i = 0; i < passband; i++)
	{
		double real = 0.0;
		double imag = -1.0 * om[i] * tau;
		DRe[i] = exp(real) * cos(imag);
		DIm[i] = exp(real) * sin(imag);
	}
	for (i = 0; i < passband; i++)
		W[i] = 1.0;
	for (i = passband; i < gridLength; i++)
		W[i] = 100.0;
	EquationErrorIIR initSolution;
	InitEquationErrorIIR(&initSolution, M, N, gridLength);
	const int iterEqnErr = 2;
	eqnerror(&initSolution, om, DRe, DIm, W, iterEqnErr);
	printf("b=[");
	for (i = 0; i < M + 1; i++)
		printf("%1.17lf ", initSolution.b[i]);
	printf("];a=[");
	for (i = 0; i < N + 1; i++)
		printf("%1.17lf ", initSolution.a[i]);
	printf("];\n");
	// Convert to SOS form
	double *sos;
	int numSecs = tf2sos(initSolution.b, M + 1, initSolution.a, N + 1, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (int j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	// Free memory
	free(sos);
	EquationErrorIIRFree(&initSolution);
	free(om);
	free(DRe);
	free(DIm);
	free(W);
	return 0;
}
// Passband linear phase bandpass filter
/*int main()
{
	int i;
	const int M = 25;
	const int N = 8;
	const double tau = 30.0;
	const int gridLength = 92;
	int passband = 10;
	double *om = (double*)malloc(gridLength * sizeof(double));
	double *DRe = (double*)malloc(gridLength * sizeof(double));
	double *DIm = (double*)malloc(gridLength * sizeof(double));
	memset(DRe, 0, gridLength * sizeof(double));
	memset(DIm, 0, gridLength * sizeof(double));
	double *W = (double*)malloc(gridLength * sizeof(double));
	linspace(0.0, 0.34, 36, om);
	linspace(0.4, 0.5, passband, om + 36);
	linspace(0.56, 1.0, 46, om + 36 + passband);
	for (i = 0; i < gridLength; i++)
		om[i] *= M_PI;
	for (i = 36; i < 36 + passband; i++)
	{
		double real = 0.0;
		double imag = -1.0 * om[i] * tau;
		DRe[i] = exp(real) * cos(imag);
		DIm[i] = exp(real) * sin(imag);
	}
	for (i = 0; i < 36; i++)
		W[i] = 10.0;
	for (i = 36; i < 36 + passband; i++)
		W[i] = 100.0;
	for (i = 36 + passband; i < gridLength; i++)
		W[i] = 10.0;
	EquationErrorIIR initSolution;
	InitEquationErrorIIR(&initSolution, M, N, gridLength);
	const int iterEqnErr = 2;
	eqnerror(&initSolution, om, DRe, DIm, W, iterEqnErr);
	printf("fvtool([");
	for (i = 0; i < M + 1; i++)
		printf("%1.17lf ", initSolution.b[i]);
	printf("],[");
	for (i = 0; i < N + 1; i++)
		printf("%1.17lf ", initSolution.a[i]);
	printf("]);\n");
	// Convert to SOS form
	double *sos;
	int numSecs = tf2sos(initSolution.b, M + 1, initSolution.a, N + 1, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (int j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	// Free memory
	free(sos);
	EquationErrorIIRFree(&initSolution);
	free(om);
	free(DRe);
	free(DIm);
	free(W);
	return 0;
}*/
// Passband linear phase bandpass & highpass filter
/*int main()
{
	int i;
	const int M = 25;
	const int N = 14;
	const int gridLength = 88;
	double *om = (double*)malloc(gridLength * sizeof(double));
	double *DRe = (double*)malloc(gridLength * sizeof(double));
	double *DIm = (double*)malloc(gridLength * sizeof(double));
	memset(DRe, 0, gridLength * sizeof(double));
	memset(DIm, 0, gridLength * sizeof(double));
	double *W = (double*)malloc(gridLength * sizeof(double));
	linspace(0.0, 0.2, 36, om);
	linspace(0.28, 0.34, 10, om + 36);
	linspace(0.43, 0.7, 26, om + 36 + 10);
	linspace(0.76, 1.0, 16, om + 36 + 10 + 26);
	for (i = 0; i < gridLength; i++)
		om[i] *= M_PI;
	double delay1 = 28.0;
	for (i = 36; i < 36 + 10; i++)
	{
		double real = 0.0;
		double imag = -om[i] * delay1;
		DRe[i] = exp(real) * cos(imag);
		DIm[i] = exp(real) * sin(imag);
	}
	double delay2 = 26.0;
	for (i = 36 + 10 + 26; i < gridLength; i++)
	{
		double real = 0.0;
		double imag = -om[i] * delay2;
		DRe[i] = exp(real) * cos(imag);
		DIm[i] = exp(real) * sin(imag);
	}
	for (i = 0; i < 36; i++)
		W[i] = 10.0;
	for (i = 36; i < 36 + 10; i++)
		W[i] = 100.0;
	for (i = 36 + 10; i < 36 + 10 + 26; i++)
		W[i] = 10.0;
	for (i = 36 + 10 + 26; i < gridLength; i++)
		W[i] = 100.0;
	EquationErrorIIR initSolution;
	InitEquationErrorIIR(&initSolution, M, N, gridLength);
	const int iterEqnErr = 2;
	eqnerror(&initSolution, om, DRe, DIm, W, iterEqnErr);
	printf("fvtool([");
	for (i = 0; i < M + 1; i++)
		printf("%1.17lf ", initSolution.b[i]);
	printf("],[");
	for (i = 0; i < N + 1; i++)
		printf("%1.17lf ", initSolution.a[i]);
	printf("]);\n");
	// Convert to SOS form
	double *sos;
	int numSecs = tf2sos(initSolution.b, M + 1, initSolution.a, N + 1, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (int j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	// Free memory
	free(sos);
	EquationErrorIIRFree(&initSolution);
	free(om);
	free(DRe);
	free(DIm);
	free(W);
	return 0;
}*/
// Group delay IIR filter (Allpass filter)
/*int main()
{
	int i;
	const int M = 8;
	const int N = 8;
	const int gridLength = 16;
	double *om = (double*)malloc(gridLength * sizeof(double));
	linspace(0.0, 1.0, gridLength, om);
	for (i = 0; i < gridLength; i++)
		om[i] *= M_PI;
	double *W = (double*)malloc(gridLength * sizeof(double));
	for (i = 0; i < gridLength; i++)
		W[i] = 1.0;
	const double phRe[16] = { 0.999999999950737,0.945743258389046,0.941157214389909,0.909809309158055,-0.227591758407713,-0.971796125574726,-0.888320467971764,-0.541219446739069,-0.162711276430479,0.171833314590005,0.444273152842578,0.655424413655085,0.811190940466005,0.917693260986536,0.979664959694248,1.0 };
	const double phIm[16] = { 0,-0.324914895286399,-0.33796907816533,0.415026530439625,-0.973756638747442,-0.235822582290091,0.459224069689757,0.840881389062365,0.986673725464686,0.985126038635165,0.89589138050509,0.755260774822106,0.584781376332966,0.397289666037134,0.20064039161462,0.0 };
	EquationErrorIIR initSolution;
	InitEquationErrorIIR(&initSolution, M, N, gridLength);
	const int iterEqnErr = 2;
	eqnerror(&initSolution, om, phRe, phIm, W, iterEqnErr);
	printf("b=[");
	for (i = 0; i < M + 1; i++)
		printf("%1.17lf ", initSolution.b[i]);
	printf("];a=[");
	for (i = 0; i < N + 1; i++)
		printf("%1.17lf ", initSolution.a[i]);
	printf("];\n");
	// Convert to SOS form
	double *sos;
	int numSecs = tf2sos(initSolution.b, M + 1, initSolution.a, N + 1, &sos);
	for (i = 0; i < numSecs; i++)
	{
		for (int j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	// Free memory
	free(sos);
	EquationErrorIIRFree(&initSolution);
	free(om);
	free(W);
	return 0;
}*/