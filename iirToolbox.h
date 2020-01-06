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
#ifndef IIRTOOLBOX_H
#define IIRTOOLBOX_H
typedef struct
{
	int num, denom, reqLen, A0xLen, twoReqLen;
	size_t memSize;
	double *memoryBuffer;
	double *b, *a;
} EquationErrorIIR;
extern void InitEquationErrorIIR(EquationErrorIIR *iir, int num, int denom, int reqLen);
extern void EquationErrorIIRFree(EquationErrorIIR *iir);
extern void eqnerror(EquationErrorIIR *iir, double *om, double *DReal, double *DImag, double *W, int iter);
extern int cplxpair(double *xRe, double *xIm, unsigned int xLen, double *sortedRe, double *sortedIm);
extern int zp2sos(double *zRe, double *zIm, unsigned int zLen, double *pRe, double *pIm, unsigned int pLen, double *sos);
extern int tf2sos(double *b, unsigned int bLen, double *a, unsigned int aLen, double **sos);
extern int designMinimumPhaseArbIIR(int gridLen, double *ff, double *aa, int gVLen, double *b, int M, double *a, int N, int iterEqnErr);
#endif