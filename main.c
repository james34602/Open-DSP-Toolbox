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
#include "utility.h"
#include "iirToolbox.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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
/*int main()
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
}*/
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
#define FREQ_VECTOR_LEN 12
int main()
{
	int i;
	double fs = 48000.0;
	double freq[FREQ_VECTOR_LEN] = { 0.0, 100.0, 400.0, 630.0, 1000.0, 1600.0, 2500.0, 4000.0, 6300.0, 10000.0, 16000.0, 24000.0 };
	double gain[FREQ_VECTOR_LEN] = { 0.0, 3.0, -2.0, 4.0, -1.0, 2.0, 4.0, 0.0, -2.0, 1.0, -0.5 };
	// To improve gain curve fitting, gain at DC is equal to first non DC bin, gain at Nyquist equal to last non Nyquist bin
	gain[0] = gain[1];
	gain[FREQ_VECTOR_LEN - 1] = gain[FREQ_VECTOR_LEN - 2];
	// Normalize to [0 ... 1]
	for (int i = 0; i < FREQ_VECTOR_LEN; i++)
		freq[i] /= (fs * 0.5);
	// Filter order
	int M = 23;
	int N = 23;
	// Allocate memory for transfer function
	double *b = (double*)malloc((M + 1) * sizeof(double));
	double *a = (double*)malloc((N + 1) * sizeof(double));
	// Filter design...designMinimumPhaseArbIIR(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
	designMinimumPhaseArbIIR(16384, freq, gain, FREQ_VECTOR_LEN, b, M, a, N, 1);
	// Print filter
	printf("fvtool([");
	for (i = 0; i < M + 1; i++)
		printf("%1.17lf ", b[i]);
	printf("],[");
	for (i = 0; i < N + 1; i++)
		printf("%1.17lf ", a[i]);
	printf("]);\n");
	// Convert to SOS form
	double *sos;
	int numSecs = tf2sos(b, M + 1, a, N + 1, &sos);
	printf("Result SOS matrix:\n");
	for (i = 0; i < numSecs; i++)
	{
		for (int j = 0; j < 6; j++)
			printf("%1.14lf ", sos[i * 6 + j]);
		printf("\n");
	}
	// Free memory
	free(sos);
	free(b);
	free(a);
	return 0;
}