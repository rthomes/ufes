/***************************************************************************
*   Copyright (C) 2005 by Fernando Cesar Meira Menandro                   *
*   fcmm@npd.ufes.br                                                      *
*																		  *
*	Created: 12-Dec-18	Renan Lima Thomes, renanlthomes@hotmail.com       *
*	Supervised by:		Fernando Cesar Meira Menandro                     *
*																		  *
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

#include "elpol_quad_2d8n.h"

elpol_quad_2D8N::elpol_quad_2D8N()
: elpol_quad_2d(nno)
{
}

elpol_quad_2D8N::~elpol_quad_2D8N()
{
}

int elpol_quad_2D8N::qnno(){
	return nno;
}

void elpol_quad_2D8N::funcao_Forma(double r, double s, double *N, double *dn){
	// RENAN
	// CALCULO DAS FUNCOES DE FORMA (N) E SUAS DERIVADAS (dn)
	// A. TABARRAEI and N. SUKUMAR, Application of polygonal finite elements in linear elasticity
	int i, j;
	double A, B, C, D, E, F, G, H, a[8], b;
	double e[14];
	const double raiz2 = sqrt(2);
	// 
	A = raiz2*(-1 + r - s);
	B = raiz2*(1 + r - s);
	C = raiz2*(-1 + r + s);
	D = raiz2*(1 + r + s);
	// Express�es que se repetem
	e[0] = (-2 * r + A);
	e[1] = (2 * s + B);
	e[2] = (-2 * r + C);
	e[3] = (-2 * s + C);
	e[4] = (-2 * r + D);
	e[5] = (-2 * s + D);
	e[8] = (2 * s + A);
	e[9] = (-2 * r + B);
	e[10] = (2 * r - A);
	e[6] = (10 * r - 40 * pow(r, 3) - 40 * r*pow(s, 2) + 14 * r*raiz2*(pow(r, 2) + pow(s, 2)) +
		2 * r*raiz2*(-3 + 7 * pow(r, 2) + 7 * pow(s, 2)));
	//e[7] = 64.*pow(b / 64, 2); Apos a vari�vel b
	e[11] = (10 * s - 40 * pow(r, 2)*s - 40 * pow(s, 3) + 14 * raiz2*s*(pow(r, 2) + pow(s, 2)) +
		2 * raiz2*s*(-3 + 7 * pow(r, 2) + 7 * pow(s, 2)));
	e[12] = (-2 + raiz2);
	e[13] = (2 - raiz2);
	// 
	a[0] = e[0]* e[1]* e[2]* e[3]* e[4]* e[5];
	a[1] = e[0]* e[8]* e[1]* e[3]* e[4]* e[5];
	a[2] = e[0]* e[9]* e[8]* e[1]* e[3]* e[4];
	a[3] = -e[0]* e[9]* e[8]* e[2]* e[3]* e[4];
	a[4] = e[9]* e[8]* e[2]* e[3]* e[4]* e[5];
	a[5] = e[9]* e[8]* e[1]* e[2]* e[3]* e[5];
	a[6] = e[0]* e[9]* e[8]* e[1]* e[2]* e[5];
	a[7] = -e[0]* e[9]* e[1]* e[2]* e[4]* e[5];
	b = 64.*(-1 + 5*pow(r, 2) - 10*pow(r, 4) + 5*pow(s, 2) - 20*pow(r, 2)*pow(s, 2) - 10*pow(s, 4) + 
		raiz2*(pow(r, 2) + pow(s, 2))*(-3 + 7*pow(r, 2) + 7*pow(s, 2)));
	e[7] = 64.*pow(b / 64, 2);
	// FUNCOES DE FORMA --
	for (i = 0; i < nno; i++){
		N[i] = a[i] / b;
	}
	// DERIVADAS ---------
	// Derivada de Ni em relacao a r
	dn[0] = -(e[0] * e[1] * e[2] * e[3] * e[4] * e[5] * e[6]) / e[7] + 	((raiz2*e[0] * e[1] * e[2] * e[3] * e[4]) + 
		(e[12] * e[0] * e[1] * e[2] * e[3] * e[5]) + (raiz2*e[0] * e[1] * e[2] * e[4] * e[5]) + 
		(e[12] * e[0] * e[1] * e[3] * e[4] * e[5]) + (raiz2*e[0] * e[2] * e[3] * e[4] * e[5]) + 
		(e[12] * e[1] * e[2] * e[3] * e[4] * e[5])) / b;
	dn[2] = -(e[0] * e[8] * e[1] * e[3] * e[4] * e[5] * e[6]) / e[7] + ((raiz2*e[0] * e[8] * e[1] * e[3] * e[4]) + 
		(e[12] * e[0] * e[8] * e[1] * e[3] * e[5]) + (raiz2*e[0] * e[8] * e[1] * e[4] * e[5]) + 
		(raiz2*e[0] * e[8] * e[3] * e[4] * e[5]) + (raiz2*e[0] * e[1] * e[3] * e[4] * e[5]) + 
		(e[12] * e[8] * e[1] * e[3] * e[4] * e[5])) / b;
	dn[4] = -(e[0] * e[9] * e[8] * e[1] * e[3] * e[4] * e[6]) / e[7] + ((e[12] * e[0] * e[9] * e[8] * e[1] * e[3]) + 
		(raiz2*e[0] * e[9] * e[8] * e[1] * e[4]) + (raiz2*e[0] * e[9] * e[8] * e[3] * e[4]) + 
		(raiz2*e[0] * e[9] * e[1] * e[3] * e[4]) + (e[12] * e[0] * e[8] * e[1] * e[3] * e[4]) + 
		(e[12] * e[9] * e[8] * e[1] * e[3] * e[4])) / b;
	dn[6] = -(e[10] * e[9] * e[8] * e[2] * e[3] * e[4] * e[6]) / e[7] + ((e[12] * e[10] * e[9] * e[8] * e[2] * e[3]) + 
		(raiz2*e[10] * e[9] * e[8] * e[2] * e[4]) + (e[12] * e[10] * e[9] * e[8] * e[3] * e[4]) + 
		(raiz2*e[10] * e[9] * e[2] * e[3] * e[4]) + (e[12] * e[10] * e[8] * e[2] * e[3] * e[4]) + 
		(e[13] * e[9] * e[8] * e[2] * e[3] * e[4])) / b;
	dn[8] = -(e[9] * e[8] * e[2] * e[3] * e[4] * e[5] * e[6]) / e[7] + (raiz2*e[9] * e[8] * e[2] * e[3] * e[4]) / b + 
		(e[12] * e[9] * e[8] * e[2] * e[3] * e[5]) / b + (raiz2*e[9] * e[8] * e[2] * e[4] * e[5]) / b + 
		(e[12] * e[9] * e[8] * e[3] * e[4] * e[5]) / b + (raiz2*e[9] * e[2] * e[3] * e[4] * e[5]) / b + 
		(e[12] * e[8] * e[2] * e[3] * e[4] * e[5]) / b;
	dn[10] = -(e[9] * e[8] * e[1] * e[2] * e[3] * e[5] * e[6]) / e[7] + (raiz2*e[9] * e[8] * e[1] * e[2] * e[3]) / b + 
		(raiz2*e[9] * e[8] * e[1] * e[2] * e[5]) / b + (e[12] * e[9] * e[8] * e[1] * e[3] * e[5]) / b + 
		(raiz2*e[9] * e[8] * e[2] * e[3] * e[5]) / b + (raiz2*e[9] * e[1] * e[2] * e[3] * e[5]) / b + 
		(e[12] * e[8] * e[1] * e[2] * e[3] * e[5]) / b;
	dn[12] = -(e[0] * e[9] * e[8] * e[1] * e[2] * e[5] * e[6]) / e[7] + (raiz2*e[0] * e[9] * e[8] * e[1] * e[2]) / b + 
		(e[12] * e[0] * e[9] * e[8] * e[1] * e[5]) / b + (raiz2*e[0] * e[9] * e[8] * e[2] * e[5]) / b + 
		(raiz2*e[0] * e[9] * e[1] * e[2] * e[5]) / b + (e[12] * e[0] * e[8] * e[1] * e[2] * e[5]) / b + 
		(e[12] * e[9] * e[8] * e[1] * e[2] * e[5]) / b;
	dn[14] = -(e[10] * e[9] * e[1] * e[2] * e[4] * e[5] * e[6]) / e[7] + (raiz2*e[10] * e[9] * e[1] * e[2] * e[4]) / b + 
		(e[12] * e[10] * e[9] * e[1] * e[2] * e[5]) / b + (e[12] * e[10] * e[9] * e[1] * e[4] * e[5]) / b + 
		(raiz2*e[10] * e[9] * e[2] * e[4] * e[5]) / b + (e[12] * e[10] * e[1] * e[2] * e[4] * e[5]) / b + 
		(e[13] * e[9] * e[1] * e[2] * e[4] * e[5]) / b;
	// Derivada de Ni em relacao a s
	dn[1] = -(e[0] * e[1] * e[2] * e[3] * e[4] * e[5] * e[11]) / e[7] + ((e[12] * e[0] * e[1] * e[2] * e[3] * e[4]) + 
		(raiz2*e[0] * e[1] * e[2] * e[3] * e[5]) + (e[12] * e[0] * e[1] * e[2] * e[4] * e[5]) + 
		(raiz2*e[0] * e[1] * e[3] * e[4] * e[5]) + (e[13] * e[0] * e[2] * e[3] * e[4] * e[5]) - 
		(raiz2*e[1] * e[2] * e[3] * e[4] * e[5])) / b;
	dn[3] = -(e[0] * e[8] * e[1] * e[3] * e[4] * e[5] * e[11]) / e[7] + ((e[12] * e[0] * e[8] * e[1] * e[3] * e[4]) + 
		(raiz2*e[0] * e[8] * e[1] * e[3] * e[5]) + (e[12] * e[0] * e[8] * e[1] * e[4] * e[5]) + 
		(e[13] * e[0] * e[8] * e[3] * e[4] * e[5]) + (e[13] * e[0] * e[1] * e[3] * e[4] * e[5]) - 
		(raiz2*e[8] * e[1] * e[3] * e[4] * e[5])) / b;
	dn[5] = -(e[0] * e[9] * e[8] * e[1] * e[3] * e[4] * e[11]) / e[7] + ((raiz2*e[0] * e[9] * e[8] * e[1] * e[3]) + 
		(e[12] * e[0] * e[9] * e[8] * e[1] * e[4]) + (e[13] * e[0] * e[9] * e[8] * e[3] * e[4]) + 
		(e[13] * e[0] * e[9] * e[1] * e[3] * e[4]) - (raiz2*e[0] * e[8] * e[1] * e[3] * e[4]) - 
		(raiz2*e[9] * e[8] * e[1] * e[3] * e[4])) / b;
	dn[7] = -(e[10] * e[9] * e[8] * e[2] * e[3] * e[4] * e[11]) / e[7] + ((raiz2*e[10] * e[9] * e[8] * e[2] * e[3]) + 
		(e[12] * e[10] * e[9] * e[8] * e[2] * e[4]) + (raiz2*e[10] * e[9] * e[8] * e[3] * e[4]) + 
		(e[13] * e[10] * e[9] * e[2] * e[3] * e[4]) - (raiz2*e[10] * e[8] * e[2] * e[3] * e[4]) +
		(raiz2*e[9] * e[8] * e[2] * e[3] * e[4])) / b;
	dn[9] = -(e[9] * e[8] * e[2] * e[3] * e[4] * e[5] * e[11]) / e[7] + ((e[12] * e[9] * e[8] * e[2] * e[3] * e[4]) + 
		(raiz2*e[9] * e[8] * e[2] * e[3] * e[5]) + (e[12] * e[9] * e[8] * e[2] * e[4] * e[5]) +
		(raiz2*e[9] * e[8] * e[3] * e[4] * e[5]) + (e[13] * e[9] * e[2] * e[3] * e[4] * e[5]) - 
		(raiz2*e[8] * e[2] * e[3] * e[4] * e[5])) / b;
	dn[11] = -(e[9] * e[8] * e[1] * e[2] * e[3] * e[5] * e[11]) / e[7] + ((e[12] * e[9] * e[8] * e[1] * e[2] * e[3]) + 
		(e[12] * e[9] * e[8] * e[1] * e[2] * e[5]) + (raiz2*e[9] * e[8] * e[1] * e[3] * e[5]) + 
		(e[13] * e[9] * e[8] * e[2] * e[3] * e[5]) + (e[13] * e[9] * e[1] * e[2] * e[3] * e[5]) - 
		(raiz2*e[8] * e[1] * e[2] * e[3] * e[5])) / b;
	dn[13] = -(e[0] * e[9] * e[8] * e[1] * e[2] * e[5] * e[11]) / e[7] + ((e[12] * e[0] * e[9] * e[8] * e[1] * e[2]) + 
		(raiz2*e[0] * e[9] * e[8] * e[1] * e[5]) + (e[13] * e[0] * e[9] * e[8] * e[2] * e[5]) + 
		(e[13] * e[0] * e[9] * e[1] * e[2] * e[5]) - (raiz2*e[0] * e[8] * e[1] * e[2] * e[5]) -
		(raiz2*e[9] * e[8] * e[1] * e[2] * e[5])) / b;
	dn[15] = -(e[10] * e[9] * e[1] * e[2] * e[4] * e[5] * e[11]) / e[7] + ((e[12] * e[10] * e[9] * e[1] * e[2] * e[4]) + 
		(raiz2*e[10] * e[9] * e[1] * e[2] * e[5]) + (raiz2*e[10] * e[9] * e[1] * e[4] * e[5]) + 
		(e[13] * e[10] * e[9] * e[2] * e[4] * e[5]) - (raiz2*e[10] * e[1] * e[2] * e[4] * e[5]) + 
		(raiz2*e[9] * e[1] * e[2] * e[4] * e[5])) / b;
}
