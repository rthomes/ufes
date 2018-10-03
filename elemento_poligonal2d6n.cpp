/***************************************************************************
*   Copyright (C) 2005 by Fernando C�sar Meira Menandro   *
*   menandro@localhost.localdomain   *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
#include "elemento_poligonal2d6n.h"

elpol2D6N::elpol2D6N()
: elpol2d(nno, ptg, ptg_tot)
{
}

elpol2D6N::~elpol2D6N()
{
}

int elpol2D6N::qptg(){
	return ptg;
}

int elpol2D6N::qnno(){
	return nno;
}

int elpol2D6N::qptg_tot() {
	return ptg_tot;
}

void elpol2D6N::pontos_de_gauss(const int p, double *r, double *s, double *w) {
	switch (p){
		// Pontos de Gauss encontrados para o elemento de 6 nos
	case 16:
		r[0] = 0.7195160831585226;
		r[1] = 0.4978928565634403;
		r[2] = 0.3482127597486885;
		r[3] = 0.7523735013471803;
		r[4] = 0.0248049647822726;
		r[5] = 0.0447155274966833;
		r[6] = 0.9209751329699979;
		r[7] = 0.0354859892823236;
		r[8] = 0.4735730939841462;
		r[9] = 0.0454572682107589;
		r[10] = 0.2338050306438317;
		r[11] = 0.2234179783782821;
		r[12] = 0.1970500443467693;
		r[13] = 0.1124927617869400;
		r[14] = 0.4030181643750864;
		r[15] = 0.2116246099606870;
		s[0] = 0.0454572682107589;
		s[1] = 0.2338050306438317;
		s[2] = 0.2234179783782821;
		s[3] = 0.1970500443467693;
		s[4] = 0.1124927617869400;
		s[5] = 0.4030181643750864;
		s[6] = 0.0354859892797277;
		s[7] = 0.9209751329789421;
		s[8] = 0.4735730939837666;
		s[9] = 0.7195160831585226;
		s[10] = 0.4978928565634403;
		s[11] = 0.3482127597486885;
		s[12] = 0.7523735013471803;
		s[13] = 0.0248049647822726;
		s[14] = 0.0447155274966833;
		s[15] = 0.2116246099606870;
		w[0] = 0.0310158499577855;
		w[1] = 0.0747959498495329;
		w[2] = -0.0079342532597081;
		w[3] = 0.0297758459676732;
		w[4] = 0.0158269989638180;
		w[5] = 0.0387697168698726;
		w[6] = 0.0104582682352531;
		w[7] = 0.0104582682344367;
		w[8] = 0.0408307010388197;
		w[9] = 0.0310158499577855;
		w[10] = 0.0747959498495329;
		w[11] = -0.0079342532597081;
		w[12] = 0.0297758459676732;
		w[13] = 0.0158269989638180;
		w[14] = 0.0387697168698726;
		w[15] = 0.0737519414707045;
		break;
	}
}

void elpol2D6N::funcao_Forma(double r, double s, double *N, double *dn){
	// RENAN
	// CALCULO DAS FUNCOES DE FORMA (N) E SUAS DERIVADAS (dn)
	// A. TABARRAEI and N. SUKUMAR, Application of polygonal finite elements in linear elasticity
	double raiz3 = sqrt(3);
	double A, B, C, D, E, F, G, H, I, c, e, den;
	// Expres~soes que se repetem
	A = (raiz3*r - s - raiz3);
	B = (2 * s + raiz3);
	C = (3 * pow(r, 2) + 6 * r - pow(s, 2) + 3);
	D = (raiz3*r + s + raiz3);
	E = (3 * pow(r, 2) - 6 * r - pow(s, 2) + 3);
	F = (4 * pow(s, 2) - 3);
	G = (raiz3*r - s + raiz3);
	H = (raiz3 - 2 * s);
	I = (raiz3*r + s - raiz3);
	c = (6 + 6 * r);
	e = (-6 + 6 * r);
	// Denominador e numerador de N, respectivamente
	den = 18 * (pow(r, 2) + pow(s, 2) - 3);
	double num[] = {A * B * C, -D * B * E, E * F, -G * H * E, I * H * C, C * F };
	// Derivadas do denominador e numerador de N, respectivamente
	double dden[] = { 36 * r, 36 * s };
	double dnum_dr[] = { c*A*B + raiz3*B*C, -(e*D*B) - raiz3*B*E, e*F,
		-(e*H*G) - raiz3*H*E, c*H*I + raiz3*H*C, c*F };
	double dnum_ds[] = { -2 * A*s*B + 2 * A*C - B*C, 2 * s*D*B - 2 * D*E - B*E, 8 * s*E - 2 * s*F,
		2 * H*G*s + H*E + 2 * G*E, -2 * H*s*I + H*C - 2 * I*C, 8 * s*C - 2 * s*F };

	// Calculo N e dn
	int j = 0;
	for (int i = 0; i < nno; i++){
		N[i] = num[i] / den;
		dn[j] = (den * dnum_dr[i] - num[i] * dden[0]) / pow(den, 2);	// Derivada de N em relacao a r
		dn[j + 1] = (den * dnum_ds[i] - num[i] * dden[1]) / pow(den, 2);	// Derivada de N em relacao a s

		j = j + 2;
	}
}
