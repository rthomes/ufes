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

void elpol2D5N::monta_n()
{
#ifdef ALEATORIO
	aleatorio
#else
	double
#endif
		// Inicializacao de variaveis ---
		r, s, J[2][2], invJ[2][2];
	double r1[2], r2[2];
	double detJ1;
	const double pi = 3.14159265358979323846;
	int i, j, n;
	i = j = n = 0;
	for (i = 0; i < 2; i++)
		for (n = 0; n < qnno(); n++)
			dn[2 * n + i] = dN[2 * n + i] = 0.0;
	// -------------------------------

	// Mudança de coordenadas dos pontos de Gauss. Coordenadas do triângulo
	// padrao para o elemento padrao ----------------------------------------
	// Pontos das extremidades do triangulo (coordenadas do elemento padrao)
	for (i = 0; i < 2; i++) {
		r1[i] = cos(2 * pi*(tri + i + 1) / qnno());
		r2[i] = sin(2 * pi*(tri + i + 1) / qnno());
	}
	detJ1 = r1[0] * r2[1] - r2[0] * r1[1];
	// Esse determinante nao muda conforme o ponto de Gauss, so depende do triangulo

	// Pontos de Gauss estão no sistema de coordenadas do triangulo,
	//entao, ele sao transformados para o sistema de coordenadas do elemento.
	// r e s: pontos de Gauss nas coordenadas do elemento.
	r = r1[0] * rpg[pg] + r1[1] * spg[pg];
	s = r2[0] * rpg[pg] + r2[1] * spg[pg];
	peso = wpg[pg];
	// -----------------------------------------------------------------------

	// Calcula N e dn para os pontos r e s (geralmente, os pontos de Gauss nas coordenadas do elemento)
	funcao_Forma(r, s, N, dn);

	// Matriz Jacobiana e Jacobiano
	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			for (n = 0; n < qnno(); n++) {
				J[i][j] += dn[2 * n + i] * this->pno[n]->qx(j);
			}
		}
	}
	detJ = J[0][0] * J[1][1] - J[1][0] * J[0][1];

	//if (detJ > 100 || detJ < 0.001){
	//	ofstream myfile;
	//	string str;
	//	str = "Entradas e Saídas/Debug_" + to_string(qnno()) + "n.txt";
	//	myfile.open(str);			
	//	myfile << "detJ = " << detJ << "\npg = " << pg << "\ntri = " << tri;
	//	myfile.close();
	//}

	invJ[0][0] = J[1][1] / detJ;
	invJ[1][1] = J[0][0] / detJ;
	invJ[0][1] = -J[0][1] / detJ;
	invJ[1][0] = -J[1][0] / detJ;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (n = 0; n < qnno(); n++)
				dN[2 * n + i] += invJ[i][j] * dn[2 * n + j];
	peso *= abs(detJ) * abs(detJ1);	// peso = peso * (detJ * detJ1)
}

void elpol2D5N::monta_rigidez()
{
#ifdef ALEATORIO
	aleatorio *xx, *yy;
	xx = new aleatorio[qnno()];
	yy = new aleatorio[qnno()];
#else
	double *xx, *yy;
	xx = new double[qnno()];
	yy = new double[qnno()];
#endif
	for (int i = 0; i<qnno()*qipn(); i++)
		for (int j = 0; j<qnno()*qipn(); j++)
			this->k[qnno()*qipn()*i + j] = 0.0;
	pontos_de_gauss(qptg(), rpg, spg, wpg);
	for (tri = 0; tri < qnno(); tri++) {	// Para cada triangulo dentro do elemento
		for (pg = 0; pg < qptg(); pg++) {
			monta_b();
			monta_c();
			for (int i = 0; i < qnno()*qipn(); i++)
				for (int j = 0; j < qnlb(); j++)
					for (int l = 0; l < qnlb(); l++)
						for (int m = 0; m < qnno()*qipn(); m++)
							k[qnno()*qipn()*i + m] +=
							b[qnno()*qipn()*j + i] * c[qnlb()*j + l] * b[qnno()*qipn()*l + m] * peso;
		}
	}
};


#ifdef ALEATORIO
void elpol2D5N::p_processa(aleatorio *xx)
{
#else
void elpol2D5N::p_processa(double *xx)
{
#endif
	pg = qptg();
	for (int i = 0; i<qnno()*qipn(); i++)
	{
		f[i] = 0.0;
		for (int n = 0; n<qnno(); n++)
			for (int j = 0; j<qipn(); j++)
				f[j] += qk(i, n*qipn() + j)*xx[qno(n)*qipn() + j];
	}
	for (int n = 0; n<qnno(); n++)
		for (int i = 0; i<qipn(); i++)
			x[n*qipn() + i] = xx[qno(n)*qipn() + i];

	// Calculo da area (Ae) e do centro (ptm) do elemento
	area_centro();
	// Zera tensao media
	for (int i = 0; i < qnlb(); i++)
		tenM[i] = 0;
	for (tri = 0; tri < qnno(); tri++) {
		//for (pg = 0; pg < lpg; pg++)
		for (pg = 0; pg < qptg(); pg++)
		{
			monta_b();
			// Calculo das coordenados dos pontos de Gauss no dominio "real"
			// e deslocamentos nos pontos de Gauss
			Peso[pg + tri * qptg()] = peso;
			ptx[pg + tri * qptg()] = pty[pg + tri * qptg()] = 0;
			for (int i = 0; i < qipn(); i++)
				des[(pg + tri * qptg())*qipn() + i] = 0;
			for (int n = 0; n < qnno(); n++) {
				ptx[pg + tri * qptg()] += N[n] * pno[n]->qx(0);
				pty[pg + tri * qptg()] += N[n] * pno[n]->qx(1);
				des[(pg + tri * qptg())*qipn()] += N[n] * x[n*qipn()];
				des[(pg + tri * qptg())*qipn() + 1] += N[n] * x[n*qipn() + 1];
			}
			//
			for (int i = 0; i < qnlb(); i++)
			{
				// Adicionei o  + tri*qptg()*qnlb() para computar os pontos de Gauss em cada triângulo
				// para 3 pontos de Gauss por triangulo, um elemento de 5 nós terá 15 pontos de Gauss.
				// Mas ainda falta escrever no arquivo de saída todos os pontos de Gauss.
				def[pg*qnlb() + i + tri * qptg()*qnlb()] = ten[pg*qnlb() + i + tri * qptg()*qnlb()] = 0;
				for (int j = 0; j < qnno()*qipn(); j++)
					def[pg*qnlb() + i + tri * qptg()*qnlb()] += b[i*qnno()*qipn() + j] * x[j];
			}
			monta_c();
			for (int i = 0; i < qnlb(); i++)
			{
				for (int j = 0; j < qnlb(); j++)
					// Revisar se é def[...+j+...]
					ten[pg*qnlb() + i + tri * qptg()*qnlb()] += c[i*qnlb() + j] * def[pg*qnlb() + j + tri * qptg()*qnlb()];
				tenM[i] += ten[pg*qnlb() + i + tri * qptg()*qnlb()] * peso;
			}
		}
	}
	for (int i = 0; i < qnlb(); i++)
		tenM[i] = tenM[i] / Ae;

	//// Tensao media
	//double lpg = ptg*qnno();
	//for (int i = 0; i < qnlb(); i++){
	//	tenM[i] = 0;
	//	for (pg = 0; pg < lpg; pg++)
	//		tenM[i] += ten[pg*qnlb() + i];
	//	tenM[i] = tenM[i] / lpg;
	//}
};