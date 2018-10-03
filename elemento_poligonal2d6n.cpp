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
		// Pontos de Gauss-Legendre encontrados para o hexagono
	case 15:
		r[0] = 0.2433678800209670;
		r[1] = -0.3638627911369946;
		r[2] = 0.1204949111160271;
		r[3] = 0.6815511544379341;
		r[4] = -0.1919704096443466;
		r[5] = -0.4895807447935853;
		r[6] = 1.0439220088499217;
		r[7] = -0.3932121491786226;
		r[8] = -0.6507098596713032;
		r[9] = 0.7376908434658901;
		r[10] = -0.6420307742629977;
		r[11] = -0.0956600692028959;
		r[12] = 0.3729917195782382;
		r[13] = -0.8321590642932697;
		r[14] = 0.4591673447150312;
		s[0] = 0.2796440497798506;
		s[1] = 0.0709407416733938;
		s[2] = -0.3505847914532461;
		s[3] = -0.1718254071120293;
		s[4] = 0.6761533172778766;
		s[5] = -0.5043279101658500;
		s[6] = -0.1486663724686596;
		s[7] = 0.9783961654680452;
		s[8] = -0.8297297929993829;
		s[9] = 0.3154472736437755;
		s[10] = 0.4811353737587404;
		s[11] = -0.7965826474025182;
		s[12] = 0.7455476498526157;
		s[13] = -0.0497535203703157;
		s[14] = -0.6957941294823040;
		w[0] = 0.3444822273297101;
		w[1] = 0.3444822273297101;
		w[2] = 0.3444822273297101;
		w[3] = 0.2327014949571080;
		w[4] = 0.2327014949571080;
		w[5] = 0.2327014949571080;
		w[6] = 0.0113025480286300;
		w[7] = 0.0113025480286300;
		w[8] = 0.0113025480286300;
		w[9] = 0.1347295079783832;
		w[10] = 0.1347295079783832;
		w[11] = 0.1347295079783832;
		w[12] = 0.1428096254906091;
		w[13] = 0.1428096254906091;
		w[14] = 0.1428096254906091;
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

void elpol2D6N::monta_n()
{
#ifdef ALEATORIO
	aleatorio
#else
	double
#endif
		// Inicializacao de variaveis ----
		r, s, J[2][2], invJ[2][2];
	int i, j, n;
	i = j = n = 0;
	for (i = 0; i < 2; i++)
		for (n = 0; n < qnno(); n++)
			dn[2 * n + i] = dN[2 * n + i] = 0.0;
	// Pontos de Gauss ---------------
	r = rpg[pg];
	s = spg[pg];
	peso = wpg[pg];

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
	peso *= abs(detJ);	// peso = peso * (detJ)
}

void elpol2D6N::monta_rigidez()
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
};


#ifdef ALEATORIO
void elpol2D5N::p_processa(aleatorio *xx)
{
#else
void elpol2D6N::p_processa(double *xx)
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
	//for (pg = 0; pg < lpg; pg++)
	for (pg = 0; pg < qptg(); pg++)
	{
		monta_b();
		// Calculo das coordenados dos pontos de Gauss no dominio "real"
		// e deslocamentos nos pontos de Gauss
		Peso[pg] = peso;
		ptx[pg] = pty[pg] = 0;
		for (int i = 0; i < qipn(); i++)
			des[pg*qipn() + i] = 0;
		for (int n = 0; n < qnno(); n++) {
			ptx[pg] += N[n] * pno[n]->qx(0);
			pty[pg] += N[n] * pno[n]->qx(1);
			des[pg*qipn()] += N[n] * x[n*qipn()];
			des[pg*qipn() + 1] += N[n] * x[n*qipn() + 1];
		}
		//
		for (int i = 0; i < qnlb(); i++)
		{
			// Adicionei o  + tri*qptg()*qnlb() para computar os pontos de Gauss em cada triângulo
			// para 3 pontos de Gauss por triangulo, um elemento de 5 nós terá 15 pontos de Gauss.
			// Mas ainda falta escrever no arquivo de saída todos os pontos de Gauss.
			def[pg*qnlb() + i] = ten[pg*qnlb() + i] = 0;
			for (int j = 0; j < qnno()*qipn(); j++)
				def[pg*qnlb() + i] += b[i*qnno()*qipn() + j] * x[j];
		}
		monta_c();
		for (int i = 0; i < qnlb(); i++)
		{
			for (int j = 0; j < qnlb(); j++)
				// Revisar se é def[...+j+...]
				ten[pg*qnlb() + i] += c[i*qnlb() + j] * def[pg*qnlb() + j];
			tenM[i] += ten[pg*qnlb() + i] * peso;
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