// Jenkins-Traub complex polynomial root finder.
#include "cpoly.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static double sr, si, tr, ti, pvr, pvi, mre;
static int nn;
static double *pr, *pi, *hr, *hi, *qpr, *qpi, *qhr, *qhi, *shr, *shi;
// COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW
static void cdivid(const double ar, const double ai, const double br, const double bi, double *cr, double *ci)
{
	double r, d;
	if (br == 0 && bi == 0)
	{
		// Division by zero, c = infinity
		*cr = DBL_MAX;
		*ci = DBL_MAX;
		return;
	}
	if (fabs(br) < fabs(bi))
	{
		r = br / bi;
		d = bi + r * br;
		*cr = (ar * r + ai) / d;
		*ci = (ai * r - ar) / d;
		return;
	}
	r = bi / br;
	d = br + r * bi;
	*cr = (ar + ai * r) / d;
	*ci = (ai - ar * r) / d;
}
// MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW
static double cmod(const double r, const double i)
{
	double ar, ai;
	ar = fabs(r);
	ai = fabs(i);
	if (ar < ai)
		return ai * sqrt(1.0 + pow((ar / ai), 2.0));
	if (ar > ai)
		return ar * sqrt(1.0 + pow((ai / ar), 2.0));
	return ar * sqrt(2.0);
}
// EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
// PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV
static void polyev(const int nn, const double sr, const double si, const double pr[], const double pi[], double qr[], double qi[], double *pvr, double *pvi)
{
	int i;
	double t;

	qr[0] = pr[0];
	qi[0] = pi[0];
	*pvr = qr[0];
	*pvi = qi[0];

	for (i = 1; i <= nn; i++)
	{
		t = (*pvr) * sr - (*pvi) * si + pr[i];
		*pvi = (*pvr) * si + (*pvi) * sr + pi[i];
		*pvr = t;
		qr[i] = *pvr;
		qi[i] = *pvi;
	}
}
// COMPUTES  T = -P(S)/H(S).
// BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
static void calct(int *bol)
{
	int n;
	double hvr, hvi;
	n = nn;
	// evaluate h(s)
	polyev(n - 1, sr, si, hr, hi, qhr, qhi, &hvr, &hvi);
	*bol = cmod(hvr, hvi) <= DBL_EPSILON * 10.0 * cmod(hr[n - 1], hi[n - 1]) ? 1 : 0;
	if (!*bol)
	{
		cdivid(-pvr, -pvi, hvr, hvi, &tr, &ti);
		return;
	}
	tr = 0;
	ti = 0;
}
// CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
// BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
static void nexth(const int bol)
{
	int j, n;
	double t1, t2;
	n = nn;
	if (!bol)
	{
		for (j = 1; j < n; j++)
		{
			t1 = qhr[j - 1];
			t2 = qhi[j - 1];
			hr[j] = tr * t1 - ti * t2 + qpr[j];
			hi[j] = tr * t2 + ti * t1 + qpi[j];
		}
		hr[0] = qpr[0];
		hi[0] = qpi[0];
		return;
	}
	// If h[s] is zero replace H with qh
	for (j = 1; j < n; j++)
	{
		hr[j] = qhr[j - 1];
		hi[j] = qhi[j - 1];
	}
	hr[0] = 0;
	hi[0] = 0;
}
// BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
// QR,QI - THE PARTIAL SUMS
// MS    -MODULUS OF THE POINT
// MP    -MODULUS OF POLYNOMIAL VALUE
static double errev(const int nn, const double qr[], const double qi[], const double ms, const double mp, const double mre)
{
	double e = cmod(qr[0], qi[0]) * mre / (DBL_EPSILON + mre);
	for (int i = 0; i <= nn; i++)
		e = e * ms + cmod(qr[i], qi[i]);
	return e * (DBL_EPSILON + mre) - mp * mre;
}
// CARRIES OUT THE THIRD STAGE ITERATION.
// L3 - LIMIT OF STEPS IN STAGE 3.
// ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
//           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
// CONV    -  .TRUE. IF ITERATION CONVERGES
static void vrshft(const int l3, double *zr, double *zi, int *conv)
{
	int b, bol;
	int i, j;
	double mp, ms, omp, relstp, r1, r2, tp;

	*conv = 0;
	b = 0;
	sr = *zr;
	si = *zi;

	// Main loop for stage three
	for (i = 1; i <= l3; i++)
	{
		// Evaluate P at S and test for convergence
		polyev(nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi);
		mp = cmod(pvr, pvi);
		ms = cmod(sr, si);
		if (mp <= 20.0 * errev(nn, qpr, qpi, ms, mp, mre))
		{
			// Polynomial value is smaller in value than a bound onthe error
			// in evaluationg P, terminate the ietartion
			*conv = 1;
			*zr = sr;
			*zi = si;
			return;
		}
		if (i != 1)
		{
			if (!(b || mp < omp || relstp >= 0.05))
			{
				// Iteration has stalled. Probably a cluster of zeros. Do 5 fixed 
				// shift steps into the cluster to force one zero to dominate
				tp = relstp;
				b = 1;
				if (relstp < DBL_EPSILON) tp = DBL_EPSILON;
				r1 = sqrt(tp);
				r2 = sr * (1 + r1) - si * r1;
				si = sr * r1 + si * (1 + r1);
				sr = r2;
				polyev(nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi);
				for (j = 1; j <= 5; j++)
				{
					calct(&bol);
					nexth(bol);
				}
				omp = DBL_MAX;
				calct(&bol);
				nexth(bol);
				calct(&bol);
				if (!bol)
				{
					relstp = cmod(tr, ti) / cmod(sr, si);
					sr += tr;
					si += ti;
				}
				continue;
			}
			// Exit if polynomial value increase significantly
			if (mp * 0.1 > omp) return;
		}
		omp = mp;

		// Calculate next iterate
		calct(&bol);
		nexth(bol);
		calct(&bol);
		if (!bol)
		{
			relstp = cmod(tr, ti) / cmod(sr, si);
			sr += tr;
			si += ti;
		}
	}
}
// COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
// INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
// APPROXIMATE ZERO IF SUCCESSFUL.
// L2 - LIMIT OF FIXED SHIFT STEPS
// ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
// CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
static void fxshft(const int l2, double *zr, double *zi, int *conv)
{
	int i, j, n;
	int test, pasd, bol;
	double otr, oti, svsr, svsi;

	n = nn;
	polyev(nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi);
	test = 1;
	pasd = 0;

	// Calculate first T = -P(S)/H(S)
	calct(&bol);

	// Main loop for second stage
	for (j = 1; j <= l2; j++)
	{
		otr = tr;
		oti = ti;

		// Compute the next H Polynomial and new t
		nexth(bol);
		calct(&bol);
		*zr = sr + tr;
		*zi = si + ti;

		// Test for convergence unless stage 3 has failed once or this
		// is the last H Polynomial
		if (!(bol || !test || j == 12))
			if (cmod(tr - otr, ti - oti) < 0.5 * cmod(*zr, *zi))
			{
				if (pasd)
				{
					// The weak convergence test has been passwed twice, start the third stage
					// Iteration, after saving the current H polynomial and shift
					for (i = 0; i < n; i++)
					{
						shr[i] = hr[i];
						shi[i] = hi[i];
					}
					svsr = sr;
					svsi = si;
					vrshft(10, zr, zi, conv);
					if (*conv)
						return;
					//The iteration failed to converge. Turn off testing and restore h,s,pv and T
					test = 0;
					for (i = 0; i < n; i++)
					{
						hr[i] = shr[i];
						hi[i] = shi[i];
					}
					sr = svsr;
					si = svsi;
					polyev(nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi);
					calct(&bol);
					continue;
				}
				pasd = 1;
			}
			else
				pasd = 0;
	}
	// Attempt an iteration with final H polynomial from second stage
	vrshft(10, zr, zi, conv);
}
// CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
// POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
static void cauchy(const int nn, double pt[], double q[], double *fn_val)
{
	int i, n;
	double x, xm, f, dx, df;
	pt[nn] = -pt[nn];
	// Compute upper estimate bound
	n = nn;
	x = exp(log(-pt[nn]) - log(pt[0])) / n;
	if (pt[n - 1] != 0)
	{
		// Newton step at the origin is better, use it
		xm = -pt[nn] / pt[n - 1];
		if (xm < x) x = xm;
	}
	// Chop the interval (0,x) until f < 0
	while (1)
	{
		xm = x * 0.1;
		f = pt[0];
		for (i = 1; i <= nn; i++)
			f = f * xm + pt[i];
		if (f <= 0)
			break;
		x = xm;
	}
	dx = x;
	// Do Newton iteration until x converges to two decimal places
	while (fabs(dx / x) > 0.005)
	{
		q[0] = pt[0];
		for (i = 1; i <= nn; i++)
			q[i] = q[i - 1] * x + pt[i];
		f = q[nn];
		df = q[0];
		for (i = 1; i < n; i++)
			df = df * x + q[i];
		dx = f / df;
		x -= dx;
	}
	*fn_val = x;
}
// RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
// THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
// INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
// BASE.
// PT - MODULUS OF COEFFICIENTS OF P
static double scale(const int nn, const double pt[])
{
	int i, pexponent;
	double hi, lo, max, min, x, sc;
	double fn_val;

	// Find largest and smallest moduli of coefficients
	hi = sqrt(DBL_MAX);
	lo = DBL_MIN / DBL_EPSILON;
	max = 0;
	min = DBL_MAX;

	for (i = 0; i <= nn; i++)
	{
		x = pt[i];
		if (x > max) max = x;
		if (x != 0 && x < min) min = x;
	}

	// Scale only if there are very large or very small components
	fn_val = 1;
	if (min >= lo && max <= hi) return fn_val;
	x = lo / min;
	if (x <= 1.0)
		sc = 1.0 / (sqrt(max)* sqrt(min));
	else
	{
		sc = x;
		if (DBL_MAX / sc > max) sc = 1.0;
	}
	pexponent = (int)(log(sc) / log(DBL_RADIX) + 0.5);
	fn_val = pow(DBL_RADIX, pexponent);
	return fn_val;
}
int cpoly(double *opr, double *opi, int degree, double *zeror, double *zeroi)
{
	int cnt1, cnt2, idnn2, i, conv, j, jj, n, nm1;
	double xx, yy, cosr, sinr, xxx, zr, zi, bnd, xni, t1, t2;
	const double RADFAC = M_PI / 180.0; // Degrees-to-radians conversion factor = pi/180
	mre = 2.0 * sqrt(2.0) * DBL_EPSILON;
	xx = sqrt(0.5);
	yy = -xx;
	cosr = cos(94.0*RADFAC);
	sinr = sin(94.0*RADFAC);
	if (degree <= 0)
		return -2;
	// Remove leading zeros
	int counter = 0;
	for (i = 0; i < degree; i++)
	{
		if (opr[i] == 0.0 && opi[i] == 0.0)
			counter++;
	}
	for (int idx = 0; idx < counter; idx++)
	{
		for (i = 0; i < degree; i++)
		{
			opr[i] = opr[i + 1];
			opi[i] = opi[i + 1];
		}
	}
	degree -= counter;
	nn = degree;

	// Algorithm fails if the leading coefficient is zero, tell upper level we may need to remove leading zero?
	if (opr[0] == 0.0 && opi[0] == 0.0)
		return -1;

	// Allocate arrays
	pr = (double*)malloc((degree + 1) * sizeof(double));
	pi = (double*)malloc((degree + 1) * sizeof(double));
	hr = (double*)malloc((degree + 1) * sizeof(double));
	hi = (double*)malloc((degree + 1) * sizeof(double));
	qpr = (double*)malloc((degree + 1) * sizeof(double));
	qpi = (double*)malloc((degree + 1) * sizeof(double));
	qhr = (double*)malloc((degree + 1) * sizeof(double));
	qhi = (double*)malloc((degree + 1) * sizeof(double));
	shr = (double*)malloc((degree + 1) * sizeof(double));
	shi = (double*)malloc((degree + 1) * sizeof(double));

	// Remove the zeros at the origin if any
	while (opr[nn] == 0 && opi[nn] == 0)
	{
		idnn2 = degree - nn;
		zeror[idnn2] = 0;
		zeroi[idnn2] = 0;
		nn--;
	}

	// Make a copy of the coefficients
	for (i = 0; i <= nn; i++)
	{
		pr[i] = opr[i];
		pi[i] = opi[i];
		shr[i] = cmod(pr[i], pi[i]);
	}

	// Scale the polynomial
	bnd = scale(nn, shr);
	if (bnd != 1)
		for (i = 0; i <= nn; i++)
		{
			pr[i] *= bnd;
			pi[i] *= bnd;
		}
	while (nn > 0)
	{
		int stop = 0;
		if (nn <= 1)
		{
			cdivid(-pr[1], -pi[1], pr[0], pi[0], &zeror[degree - 1], &zeroi[degree - 1]);
			// Deallocate arrays
			free(pr);
			free(pi);
			free(hr);
			free(hi);
			free(qpr);
			free(qpi);
			free(qhr);
			free(qhi);
			free(shr);
			free(shi);
			return degree;
		}

		// Calculate bnd, alower bound on the modulus of the zeros
		for (i = 0; i <= nn; i++)
			shr[i] = cmod(pr[i], pi[i]);

		cauchy(nn, shr, shi, &bnd);

		// Outer loop to control 2 Major passes with different sequences of shifts
		for (cnt1 = 1; (cnt1 <= 2) && !stop; cnt1++)
		{
			// First stage calculation, no shift
			// COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
			// POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
			const int l1 = 5;
			n = nn;
			nm1 = n - 1;
			for (i = 0; i < n; i++)
			{
				xni = nn - i;
				hr[i] = xni * pr[i] / n;
				hi[i] = xni * pi[i] / n;
			}
			for (jj = 1; jj <= l1; jj++)
			{
				if (cmod(hr[n - 1], hi[n - 1]) > DBL_EPSILON * 10.0 * cmod(pr[n - 1], pi[n - 1]))
				{
					cdivid(-pr[nn], -pi[nn], hr[n - 1], hi[n - 1], &tr, &ti);
					for (i = 0; i < nm1; i++)
					{
						j = nn - i - 1;
						t1 = hr[j - 1];
						t2 = hi[j - 1];
						hr[j] = tr * t1 - ti * t2 + pr[j];
						hi[j] = tr * t2 + ti * t1 + pi[j];
					}
					hr[0] = pr[0];
					hi[0] = pi[0];
				}
				else
				{
					// If the constant term is essentially zero, shift H coefficients
					for (i = 0; i < nm1; i++)
					{
						j = nn - i - 1;
						hr[j] = hr[j - 1];
						hi[j] = hi[j - 1];
					}
					hr[0] = 0.0;
					hi[0] = 0.0;
				}
			}
			// Inner loop to select a shift
			for (cnt2 = 1; (cnt2 <= 9) && !stop; cnt2++)
			{
				// Shift is chosen with modulus bnd and amplitude rotated by 94 degree from the previous shif
				xxx = cosr * xx - sinr * yy;
				yy = sinr * xx + cosr * yy;
				xx = xxx;
				sr = bnd * xx;
				si = bnd * yy;

				// Second stage calculation, fixed shift
				fxshft(10 * cnt2, &zr, &zi, &conv);
				if (conv)
				{
					// The second stage jumps directly to the third stage ieration
					// If successful the zero is stored and the polynomial deflated
					idnn2 = degree - nn;
					zeror[idnn2] = zr;
					zeroi[idnn2] = zi;
					nn--;
					for (i = 0; i <= nn; i++)
					{
						pr[i] = qpr[i];
						pi[i] = qpi[i];
					}
					stop = 1;
					break;
				}
				// If the iteration is unsuccessful another shift is chosen
			}
			// if 9 shifts fail, the outer loop is repeated with another sequence of shifts
			if (stop)
				break;
		}
		// The zerofinder has failed on two major passes
		// return empty handed with the number of roots found (less than the original degree)
		if (!stop)
		{
			degree -= nn;
			// Deallocate arrays
			free(pr);
			free(pi);
			free(hr);
			free(hi);
			free(qpr);
			free(qpi);
			free(qhr);
			free(qhi);
			free(shr);
			free(shi);
			return degree;
		}
	}
	// Deallocate arrays
	free(pr);
	free(pi);
	free(hr);
	free(hi);
	free(qpr);
	free(qpi);
	free(qhr);
	free(qhi);
	free(shr);
	free(shi);
	return degree;
}