﻿/***************************************************************************
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

#include "elpol_tri_2d.h"

elpol_tri_2d::elpol_tri_2d(int nno)
: isop2d(nno, ptg, 2) // 2 indica mapeamento em 2 niveis
{
	rpg = new double[qptg()];
	spg = new double[qptg()];
	wpg = new double[qptg()];
}

elpol_tri_2d::~elpol_tri_2d()
{
}

int elpol_tri_2d::qptg(){
	return ptg;
}

void elpol_tri_2d::pontos_de_gauss(const int p, double *r, double *s, double *w) {
	switch (p){
	case 1:
		r[0] = 1.0 / 3;
		s[0] = 1.0 / 3;
		w[0] = 0.5;
		break;
	case 3:
		r[0] = 1.0 / 6;
		s[0] = 1.0 / 6;
		r[1] = 2.0 / 3;
		s[1] = 1.0 / 6;
		r[2] = 1.0 / 6;
		s[2] = 2.0 / 3;
		w[0] = 1.0 / 6;
		w[1] = 1.0 / 6;
		w[2] = 1.0 / 6;
		break;
	/*case 3:
		r[0] = 0;
		s[0] = 1.0 / 2;
		r[1] = 1.0 / 2;
		s[1] = 0;
		r[2] = 1.0 / 2;
		s[2] = 1.0 / 2;
		w[0] = 1.0 / 6;
		w[1] = 1.0 / 6;
		w[2] = 1.0 / 6;
		break;*/
	case 4:
		r[0] = 1.0 / 3;
		s[0] = 1.0 / 3;
		r[1] = 1.0 / 5;
		s[1] = 1.0 / 5;
		r[2] = 1.0 / 5;
		s[2] = 3.0 / 5;
		r[3] = 3.0 / 5;
		s[3] = 1.0 / 5;
		w[0] = -27.0 / 96;
		w[1] = 25.0 / 96;
		w[2] = 25.0 / 96;
		w[3] = 25.0 / 96;
		break;
	/*case 4:
		r[0] = 1.0 / 3;
		s[0] = 1.0 / 3;
		r[1] = 2.0 / 15;
		s[1] = 11. / 15;
		r[2] = 2.0 / 15;
		s[2] = 2.0 / 15;
		r[3] = 11. / 15;
		s[3] = 2.0 / 15;
		w[0] = -27.0 / 96;
		w[1] = 25.0 / 96;
		w[2] = 25.0 / 96;
		w[3] = 25.0 / 96;
		break;*/
	case 6:
		r[0] = 0.4459484909159650;
		r[1] = 0.4459484909159650;
		r[2] = 0.1081030181680700;
		r[3] = 0.0915762135097710;
		r[4] = 0.0915762135097710;
		r[5] = 0.8168475729804590;
		s[0] = 0.4459484909159650;
		s[1] = 0.1081030181680700;
		s[2] = 0.4459484909159650;
		s[3] = 0.0915762135097710;
		s[4] = 0.8168475729804590;
		s[5] = 0.0915762135097710;
		w[0] = 0.1116907948390055;
		w[1] = 0.1116907948390055;
		w[2] = 0.1116907948390055;
		w[3] = 0.0549758718276610;
		w[4] = 0.0549758718276610;
		w[5] = 0.0549758718276610;
		break;
	case 7:
		r[0] = 0.3333333333333330;
		r[1] = 0.4701420641051150;
		r[2] = 0.4701420641051150;
		r[3] = 0.0597158717897700;
		r[4] = 0.1012865073234560;
		r[5] = 0.1012865073234560;
		r[6] = 0.7974269853530870;
		s[0] = 0.3333333333333330;
		s[1] = 0.4701420641051150;
		s[2] = 0.0597158717897700;
		s[3] = 0.4701420641051150;
		s[4] = 0.1012865073234560;
		s[5] = 0.7974269853530870;
		s[6] = 0.1012865073234560;
		w[0] = 0.1125000000000000;
		w[1] = 0.0661970763942530;
		w[2] = 0.0661970763942530;
		w[3] = 0.0661970763942530;
		w[4] = 0.0629695902724135;
		w[5] = 0.0629695902724135;
		w[6] = 0.0629695902724135;
		break;
	case 12:
		r[0] = 0.2492867451709100;
		r[1] = 0.2492867451709100;
		r[2] = 0.5014265096581790;
		r[3] = 0.0630890144915020;
		r[4] = 0.0630890144915020;
		r[5] = 0.8738219710169960;
		r[6] = 0.6365024991213990;
		r[7] = 0.6365024991213990;
		r[8] = 0.3103524510337840;
		r[9] = 0.3103524510337840;
		r[10] = 0.0531450498448170;
		r[11] = 0.0531450498448170;
		s[0] = 0.2492867451709100;
		s[1] = 0.5014265096581790;
		s[2] = 0.2492867451709100;
		s[3] = 0.0630890144915020;
		s[4] = 0.8738219710169960;
		s[5] = 0.0630890144915020;
		s[6] = 0.3103524510337840;
		s[7] = 0.0531450498448170;
		s[8] = 0.6365024991213990;
		s[9] = 0.0531450498448170;
		s[10] = 0.6365024991213990;
		s[11] = 0.3103524510337840;
		w[0] = 0.0583931378631895;
		w[1] = 0.0583931378631895;
		w[2] = 0.0583931378631895;
		w[3] = 0.0254224531851035;
		w[4] = 0.0254224531851035;
		w[5] = 0.0254224531851035;
		w[6] = 0.0414255378091870;
		w[7] = 0.0414255378091870;
		w[8] = 0.0414255378091870;
		w[9] = 0.0414255378091870;
		w[10] = 0.0414255378091870;
		w[11] = 0.0414255378091870;
		break;
	case 16:
		r[0] = 0.333333333;
		r[1] = 0.459292588;
		r[2] = 0.459292588;
		r[3] = 0.081414823;
		r[4] = 0.170569308;
		r[5] = 0.170569308;
		r[6] = 0.658861384;
		r[7] = 0.050547228;
		r[8] = 0.050547228;
		r[9] = 0.898905543;
		r[10] = 0.26311283;
		r[11] = 0.728492393;
		r[12] = 0.008394777;
		r[13] = 0.728492393;
		r[14] = 0.26311283;
		r[15] = 0.008394777;
		s[0] = 0.333333333;
		s[1] = 0.459292588;
		s[2] = 0.081414823;
		s[3] = 0.459292588;
		s[4] = 0.170569308;
		s[5] = 0.658861384;
		s[6] = 0.170569308;
		s[7] = 0.050547228;
		s[8] = 0.898905543;
		s[9] = 0.050547228;
		s[10] = 0.728492393;
		s[11] = 0.008394777;
		s[12] = 0.26311283;
		s[13] = 0.26311283;
		s[14] = 0.008394777;
		s[15] = 0.728492393;
		w[0] = 0.144315608 / 2;
		w[1] = 0.095091634 / 2;
		w[2] = 0.095091634 / 2;
		w[3] = 0.095091634 / 2;
		w[4] = 0.103217371 / 2;
		w[5] = 0.103217371 / 2;
		w[6] = 0.103217371 / 2;
		w[7] = 0.032458498 / 2;
		w[8] = 0.032458498 / 2;
		w[9] = 0.032458498 / 2;
		w[10] = 0.027230314 / 2;
		w[11] = 0.027230314 / 2;
		w[12] = 0.027230314 / 2;
		w[13] = 0.027230314 / 2;
		w[14] = 0.027230314 / 2;
		w[15] = 0.027230314 / 2;
		break;
	case 25:
		r[0] = 0.333333333;
		r[1] = 0.485577633;
		r[2] = 0.485577633;
		r[3] = 0.028844733;
		r[4] = 0.109481575;
		r[5] = 0.109481575;
		r[6] = 0.781036849;
		r[7] = 0.550352942;
		r[8] = 0.550352942;
		r[9] = 0.307939839;
		r[10] = 0.307939839;
		r[11] = 0.141707219;
		r[12] = 0.141707219;
		r[13] = 0.728323905;
		r[14] = 0.728323905;
		r[15] = 0.246672561;
		r[16] = 0.246672561;
		r[17] = 0.025003535;
		r[18] = 0.025003535;
		r[19] = 0.923655934;
		r[20] = 0.923655934;
		r[21] = 0.066803251;
		r[22] = 0.066803251;
		r[23] = 0.009540815;
		r[24] = 0.009540815;
		s[0] = 0.333333333;
		s[1] = 0.485577633;
		s[2] = 0.028844733;
		s[3] = 0.485577633;
		s[4] = 0.109481575;
		s[5] = 0.781036849;
		s[6] = 0.109481575;
		s[7] = 0.307939839;
		s[8] = 0.141707219;
		s[9] = 0.550352942;
		s[10] = 0.141707219;
		s[11] = 0.550352942;
		s[12] = 0.307939839;
		s[13] = 0.246672561;
		s[14] = 0.025003535;
		s[15] = 0.728323905;
		s[16] = 0.025003535;
		s[17] = 0.728323905;
		s[18] = 0.246672561;
		s[19] = 0.066803251;
		s[20] = 0.009540815;
		s[21] = 0.923655934;
		s[22] = 0.009540815;
		s[23] = 0.923655934;
		s[24] = 0.066803251;
		w[0] = 0.045408995;
		w[1] = 0.0183629788782335;
		w[2] = 0.0183629788782335;
		w[3] = 0.0183629788782335;
		w[4] = 0.0226605295;
		w[5] = 0.0226605295;
		w[6] = 0.0226605295;
		w[7] = 0.0363789585;
		w[8] = 0.0363789585;
		w[9] = 0.0363789585;
		w[10] = 0.0363789585;
		w[11] = 0.0363789585;
		w[12] = 0.0363789585;
		w[13] = 0.0141636215;
		w[14] = 0.0141636215;
		w[15] = 0.0141636215;
		w[16] = 0.0141636215;
		w[17] = 0.0141636215;
		w[18] = 0.0141636215;
		w[19] = 0.0047108335;
		w[20] = 0.0047108335;
		w[21] = 0.0047108335;
		w[22] = 0.0047108335;
		w[23] = 0.0047108335;
		w[24] = 0.0047108335;
		break;
	}
}

void elpol_tri_2d::monta_n()
{
#ifdef ALEATORIO
	aleatorio
#else
	double
#endif
	// Inicializacao de variaveis ---
		r, s, J[2][2], invJ[2][2];
	double detJ1;
	int i, j, n;
	i = j = n = 0;
	// Zerar dn e dN
	for (i = 0; i < 2; i++)
	for (n = 0; n < qnno(); n++)
		dn[2 * n + i] = dN[2 * n + i] = 0.0;
	
	// Mudança de coordenadas dos pontos de Gauss. Coordenadas do triangulo
	// padrao para o elemento padrao. Atencao: peso alterado nessa funcao
	map2(&r, &s, &detJ1);

	// Calcula N e dn para os pontos r e s (geralmente, os pontos de Gauss nas coordenadas do elemento)
	funcao_Forma(r, s, N, dn);

	// Matriz Jacobiana e Jacobiano
	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			for (n = 0; n < qnno(); n++){
				J[i][j] += dn[2 * n + i] * this->pno[n]->qx(j);
			}
		}
	}
	detJ = J[0][0] * J[1][1] - J[1][0] * J[0][1];

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

void elpol_tri_2d::map2(double *r, double *s, double *detJ1) {
	// Mudança de coordenadas dos pontos de Gauss. Coordenadas do triângulo
	// padrao para o elemento padrao.
	// Poderia-se melhorar o tempo calculando todos os vertices das 
	// sub-divisoes uma unica vez
	double r1[2], r2[2];
	const double pi = 3.14159265358979323846;
	int i;
	// Pontos das extremidades do triangulo (coordenadas do elemento padrao)
	for (i = 0; i < 2; i++) {
		r1[i] = cos(2 * pi*(tri + i + 1) / qnno());
		r2[i] = sin(2 * pi*(tri + i + 1) / qnno());
	}
	*detJ1 = r1[0] * r2[1] - r2[0] * r1[1];
	// Esse determinante nao muda conforme o ponto de Gauss, so depende do triangulo

	// Pontos de Gauss estão no sistema de coordenadas do triangulo,
	// entao, ele sao transformados para o sistema de coordenadas do elemento.
	// r e s: pontos de Gauss nas coordenadas do elemento.
	*r = r1[0] * rpg[pg] + r1[1] * spg[pg];
	*s = r2[0] * rpg[pg] + r2[1] * spg[pg];
	peso = wpg[pg];
}

void elpol_tri_2d::monta_rigidez()
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
	for (tri = 0; tri < qnno(); tri++){	// Para cada triangulo dentro do elemento
		for (pg = 0; pg < qptg(); pg++){
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
void elpol_tri_2d::p_processa(aleatorio *xx)
{
#else
void elpol_tri_2d::p_processa(double *xx)
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
	for (tri = 0; tri < qnno(); tri++){
		for (pg = 0; pg < qptg(); pg++)
		{
			monta_b();
			// Calculo das coordenados dos pontos de Gauss no dominio "real"
			// e deslocamentos nos pontos de Gauss
			Peso[pg + tri*qptg()] = peso;
			ptx[pg + tri*qptg()] = pty[pg + tri*qptg()] = 0;
			for (int i = 0; i < qipn(); i++)
				des[(pg + tri*qptg())*qipn() + i] = 0;
			for (int n = 0; n < qnno(); n++){
				ptx[pg + tri*qptg()] += N[n] * pno[n]->qx(0);
				pty[pg + tri*qptg()] += N[n] * pno[n]->qx(1);
				des[(pg + tri*qptg())*qipn()] += N[n] * x[n*qipn()];
				des[(pg + tri*qptg())*qipn() + 1] += N[n] * x[n*qipn() + 1];
			}
			//
			for (int i = 0; i < qnlb(); i++)
			{
				// Adicionei o  + tri*qptg()*qnlb() para computar os pontos de Gauss em cada triângulo
				// para 3 pontos de Gauss por triangulo, um elemento de 5 nós terá 15 pontos de Gauss.
				// Mas ainda falta escrever no arquivo de saída todos os pontos de Gauss.
				def[pg*qnlb() + i + tri*qptg()*qnlb()] = ten[pg*qnlb() + i + tri*qptg()*qnlb()] = 0;
				for (int j = 0; j < qnno()*qipn(); j++)
					def[pg*qnlb() + i + tri*qptg()*qnlb()] += b[i*qnno()*qipn() + j] * x[j];
			}
			monta_c();
			for (int i = 0; i < qnlb(); i++)
			{
				for (int j = 0; j < qnlb(); j++)
					// Revisar se é def[...+j+...]
					ten[pg*qnlb() + i + tri*qptg()*qnlb()] += c[i*qnlb() + j] * def[pg*qnlb() + j + tri*qptg()*qnlb()];
				tenM[i] += ten[pg*qnlb() + i + tri*qptg()*qnlb()] * peso;
			}
		}
	}
	for (int i = 0; i < qnlb(); i++)
		tenM[i] = tenM[i] / Ae;
};