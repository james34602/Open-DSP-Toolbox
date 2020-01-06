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