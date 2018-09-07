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
#include "elemento_poligonal2d5n.h"

elpol2D5N::elpol2D5N()
: elpol2d(nno, ptg)
{
}

elpol2D5N::~elpol2D5N()
{
}

int elpol2D5N::qptg(){
	return ptg;
}

int elpol2D5N::qnno(){
	return nno;
}

void elpol2D5N::pontos_de_gauss(const int p, double *r, double *s, double *w) {
	switch (p){
		// Pontos de Gauss encontrados para o elemento de 5 nos
	case 16:
		r[0] = 0.0760673942478839;
		r[1] = 0.3247270159267784;
		r[2] = 0.8943643339137815;
		r[3] = 0.1104341136719281;
		r[4] = 0.7597491311679048;
		r[5] = 0.6460735503492643;
		r[6] = 0.6257963489099985;
		r[7] = 0.4651002824242623;
		r[8] = 0.3207916796031768;
		r[9] = 0.3275324375145142;
		r[10] = 0.0866851930002782;
		r[11] = 0.0422352459237475;
		r[12] = 0.0341036371449512;
		r[13] = 0.2228307870144855;
		r[14] = 0.0291258789079805;
		r[15] = 0.1816726588421845;
		s[0] = 0.3275324375145142;
		s[1] = 0.0866851930002782;
		s[2] = 0.0422352459237475;
		s[3] = 0.0341036371449512;
		s[4] = 0.2228307870144855;
		s[5] = 0.0291258789079805;
		s[6] = 0.1816726588421845;
		s[7] = 0.4651002824286803;
		s[8] = 0.3207916795980526;
		s[9] = 0.0760673942478839;
		s[10] = 0.3247270159267784;
		s[11] = 0.8943643339137815;
		s[12] = 0.1104341136719281;
		s[13] = 0.7597491311679048;
		s[14] = 0.6460735503492643;
		s[15] = 0.6257963489099985;
		w[0] = 0.0282442583132088;
		w[1] = 0.0284586587203032;
		w[2] = 0.0173628081849766;
		w[3] = 0.0162037480400389;
		w[4] = 0.0160852820390569;
		w[5] = 0.0263938091488609;
		w[6] = 0.0551745907121489;
		w[7] = 0.0481578165251256;
		w[8] = 0.0759959473345565;
		w[9] = 0.0282442583132088;
		w[10] = 0.0284586587203032;
		w[11] = 0.0173628081849766;
		w[12] = 0.0162037480400389;
		w[13] = 0.0160852820390569;
		w[14] = 0.0263938091488609;
		w[15] = 0.0551745907121489;
		break;
	}
}

void elpol2D5N::funcao_Forma(double r, double s, double *N, double *dn){
	// RENAN
	// Calcula N e dn para os pontos r e s (geralmente, os pontos de Gauss nas coordenadas do elemento)
	// Coeficientes para as fun��es de forma(?)
	// Expressões que se repetem
	double A, B, C, D, E, den, den2;
	A = (0.8090169943749475 + 1.*r);
	B = (-1. + 1.*r - 0.7265425280053609*s);
	C = (2.6180339887498953 + 1.*r + 3.077683537175254*s);
	D = (-1. + 1.*r + 0.7265425280053609*s);
	E = (2.6180339887498953 + 1.*r - 3.077683537175254*s);
	den = (122.60990336999416 - 17.888543819998326*pow(r, 2) - 17.88854381999832*pow(s, 2));
	den2 = pow(122.60990336999416 - 17.888543819998326*pow(r, 2) - 17.88854381999832*pow(s, 2), 2);
	// Funcao de forma
	N[0] = (-11.577708763999665*A*B*C) / den;
	N[1] = (9.366563145999496*B*D*C) / den;
	N[2] = (9.366563145999496*E*B*D) / den;
	N[3] = (-11.577708763999665*A*E*D) / den;
	N[4] = (4.422291236000336*A*E*C) / den;
	// Derivada em relação a r
	dn[0] = (-414.2167011199733*r*A*B*C) / den2 - (11.577708763999665*A*B) / den - (11.577708763999665*A*C) / den - (11.577708763999665*B*C) / den;
	dn[2] = (335.1083505599867*r*B*D*C) / den2 + (9.366563145999496*B*D) / den + (9.366563145999496*B*C) / den + (9.366563145999496*D*C) / den;
	dn[4] = (335.1083505599867*r*E*B*D) / den2 + (9.366563145999496*E*B) / den + (9.366563145999496*E*D) / den + (9.366563145999496*B*D) / den;
	dn[6] = (-414.2167011199733*r*A*E*D) / den2 - (11.577708763999665*A*E) / den - (11.577708763999665*A*D) / den - (11.577708763999665*E*D) / den;
	dn[8] = (158.21670111997315*r*A*E*C) / den2 + (4.422291236000336*A*E) / den + (4.422291236000336*A*C) / den + (4.422291236000336*E*C) / den;
	// Derivada em relação a s
	dn[1] = (-414.21670111997315*A*B*s*C) / den2 - (35.632523661171426*A*B) / den + (8.41169779390614*A*C) / den;
	dn[3] = (335.1083505599866*B*D*s*C) / den2 + (28.827317194355103*B*D) / den + (6.80520646681632*B*C) / den - (6.80520646681632*D*C) / den;
	dn[5] = (335.1083505599866*E*B*D*s) / den2 + (6.80520646681632*E*B) / den - (6.80520646681632*E*D) / den - (28.827317194355103*B*D) / den;
	dn[7] = (-414.21670111997315*A*E*D*s) / den2 - (8.41169779390614*A*E) / den + (35.632523661171426*A*D) / den;
	dn[9] = (158.21670111997307*A*E*s*C) / den2 + (13.61041293363264*A*E) / den - (13.61041293363264*A*C) / den;
}