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
void linspace(double a, double b, const int n, double *y)
{
	int i;
	double d = (b - a) / (double)(n - 1);
	for (i = 0; i < n; i++)
		y[i] = a + i * d;
}