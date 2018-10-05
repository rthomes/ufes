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
: elpol2d(nno, ptg, ptg_tot)
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

int elpol2D5N::qptg_tot() {
	return ptg_tot;
}

void elpol2D5N::pontos_de_gauss(const int p, double *r, double *s, double *w) {
	// Pontos de Gauss-Legendre encontrados para pentagono
	switch (p){
	case 16:
		r[0] = 0.5912188450226953;
		r[1] = -0.2421899246491391;
		r[2] = -0.7409004501886392;
		r[3] = -0.2157117358475382;
		r[4] = 0.6075832656626214;
		r[5] = 0.5148261490806013;
		r[6] = 0.1474010604237443;
		r[7] = -0.4237272837609504;
		r[8] = -0.4092789237486831;
		r[9] = 0.1707789980052876;
		r[10] = 0.8677983023836097;
		r[11] = 0.2624214150678483;
		r[12] = -0.7056129484958356;
		r[13] = -0.6985142001403036;
		r[14] = 0.2739074311846809;
		r[15] = 0.0000000000000003;
		s[0] = 0.4467522043916261;
		s[1] = 0.7003365585468113;
		s[2] = -0.0139204076455661;
		s[3] = -0.7089398436090251;
		s[4] = -0.4242285116838463;
		s[5] = 0.0122905091237966;
		s[6] = 0.4934267400310196;
		s[7] = 0.2926639871734320;
		s[8] = -0.3125504486747753;
		s[9] = -0.4858307876534727;
		s[10] = 0.0060385560269210;
		s[11] = 0.8271912467456078;
		s[12] = 0.5051937496582656;
		s[13] = -0.5149643385527941;
		s[14] = -0.8234592138780004;
		s[15] = 0.0000000000000001;
		w[0] = 0.1126522280160471;
		w[1] = 0.1126522280160471;
		w[2] = 0.1126522280160471;
		w[3] = 0.1126522280160471;
		w[4] = 0.1126522280160471;
		w[5] = 0.2289069524904512;
		w[6] = 0.2289069524904512;
		w[7] = 0.2289069524904512;
		w[8] = 0.2289069524904512;
		w[9] = 0.2289069524904512;
		w[10] = 0.0802467001918499;
		w[11] = 0.0802467001918499;
		w[12] = 0.0802467001918499;
		w[13] = 0.0802467001918499;
		w[14] = 0.0802467001918499;
		w[15] = 0.2686118872461438;
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