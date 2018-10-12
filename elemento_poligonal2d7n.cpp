﻿/***************************************************************************
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
#include "elemento_poligonal2d7n.h"

elpol2D7N::elpol2D7N()
: elpol2d(nno, ptg, ptg_tot)
{
}

elpol2D7N::~elpol2D7N()
{
}

int elpol2D7N::qnno(){
	return nno;
}

int elpol2D7N::qptg(){
	return ptg;
}

int elpol2D7N::qptg_tot() {
	return ptg_tot;
}

void elpol2D7N::pontos_de_gauss(const int p, double *r, double *s, double *w) {
	// Pontos de Gauss-Legendre encontrados para heptagono
	switch (p) {
		r[0] = 0.6255681659215563;
		r[1] = -0.2142726536339368;
		r[2] = -0.6255681659215564;
		r[3] = 0.2142726536339366;
		s[0] = 0.2142726536339368;
		s[1] = 0.6255681659215563;
		s[2] = -0.2142726536339366;
		s[3] = -0.6255681659215564;
		w[0] = 0.6841025471595260;
		w[1] = 0.6841025471595260;
		w[2] = 0.6841025471595260;
		w[3] = 0.6841025471595260;
		break;
	case 16:
		r[0] = 0.1276168301307205;
		r[1] = 0.1276168301307205;
		r[2] = -0.4053425121477094;
		r[3] = -0.4053425121477094;
		r[4] = -0.7879519226464595;
		r[5] = -0.7879519226464596;
		r[6] = 0.0828186085296878;
		r[7] = 0.0828186085296876;
		r[8] = 0.6127070734564268;
		r[9] = 0.6127070734564268;
		r[10] = -0.3755322104713331;
		r[11] = -0.3755322104713333;
		r[12] = 0.1539942952564272;
		r[13] = 0.1539942952564274;
		r[14] = 0.5418216463061164;
		r[15] = 0.8866988426177320;
		s[0] = -0.1490668093012800;
		s[1] = 0.1490668093012800;
		s[2] = 0.2259281525316539;
		s[3] = -0.2259281525316538;
		s[4] = 0.3095516214853822;
		s[5] = -0.3095516214853821;
		s[6] = 0.5100806420522777;
		s[7] = -0.5100806420522777;
		s[8] = 0.5510981065938946;
		s[9] = -0.5510981065938946;
		s[10] = 0.7417535573906761;
		s[11] = -0.7417535573906759;
		s[12] = 0.8687581555362361;
		s[13] = -0.8687581555362361;
		s[14] = 0.0000000000102235;
		s[15] = -0.0000000000040105;
		w[0] = 0.1185983078311412;
		w[1] = 0.1185983078311412;
		w[2] = 0.2453764631968670;
		w[3] = 0.2453764631968670;
		w[4] = 0.1447750410527415;
		w[5] = 0.1447750410527415;
		w[6] = 0.2670118477615240;
		w[7] = 0.2670118477615240;
		w[8] = 0.1681763071628871;
		w[9] = 0.1681763071628871;
		w[10] = 0.1524140290242357;
		w[11] = 0.1524140290242357;
		w[12] = 0.0870315191010727;
		w[13] = 0.0870315191010727;
		w[14] = 0.2713662505935009;
		w[15] = 0.0982769077090027;
		break;
	}
}

void elpol2D7N::funcao_Forma(double r, double s, double *N, double *dn){
	// RENAN
	// CALCULO DAS FUNCOES DE FORMA (N) E SUAS DERIVADAS (dn)
	double b, b2, dbdr, dbds, A, B, C, D, E, F, G;
	// Expressoes que se repetem	
	A = (0.9009688679024191 + 1.*r);
	B = (-4.048917339522304 + 1.*r - 4.381286267534823*s);
	C = (1.445041867912629 + 1.*r - 1.253960337662704*s);
	D = (-1. + 1.*r - 0.4815746188075287*s);
	E = (1.445041867912629 + 1.*r + 1.253960337662704*s);
	F = (-1. + 1.*r + 0.4815746188075287*s);
	G = (-4.048917339522304 + 1.*r + 4.381286267534823*s);
	b = (1412.7983817039994 + 41.27066882957638*pow(r, 4) + pow(s, 2)*(-762.7594506989002 + 41.27066882957639*pow(s, 2)) + 
		pow(r, 2)*(-762.7594506988999 + 82.54133765915262*pow(s, 2)));
	b2 = pow(b, 2);
	dbdr = (165.08267531830552*pow(r, 3) + 2 * r*(-762.7594506988999 + 82.54133765915262*pow(s, 2)));
	dbds = (165.08267531830523*pow(r, 2)*s + 82.54133765915277*pow(s, 3) + 2 * s*(-762.7594506989002 + 41.27066882957639*pow(s, 2)));
	// Funcoes de Forma
	N[0] = (26.495528883381816*A*B*C*D*E) / b;
	N[1] = (-38.28714854897507*A*B*D*F*E) / b;
	N[2] = (8.519692053642082*B*D*F*E*G) / b;
	N[3] = (8.519692053642082*B*C*D*F*G) / b;
	N[4] = (-38.28714854897507*A*C*D*F*G) / b;
	N[5] = (26.495528883381816*A*C*F*E*G) / b;
	N[6] = (6.543855223902344*A*B*C*E*G) / b;
	// Derivadas em relacao a r 
	dn[0] = (-26.495528883381816*A*B*C*D*E*dbdr) / b2 + (26.495528883381816*A*B*C*D) / b + (26.495528883381816*A*B*C*E) / b + (26.495528883381816*A*B*D*E) / b + (26.495528883381816*A*C*D*E) / b + (26.495528883381816*B*C*D*E) / b;
	dn[2] = (38.28714854897507*A*B*D*F*E*dbdr) / b2 - (38.28714854897507*A*B*D*F) / b - (38.28714854897507*A*B*D*E) / b - (38.28714854897507*A*B*F*E) / b - (38.28714854897507*A*D*F*E) / b - (38.28714854897507*B*D*F*E) / b;
	dn[4] = (-8.519692053642082*B*D*F*E*G*dbdr) / b2 + (8.519692053642082*B*D*F*E) / b + (8.519692053642082*B*D*F*G) / b + (8.519692053642082*B*D*E*G) / b + (8.519692053642082*B*F*E*G) / b + (8.519692053642082*D*F*E*G) / b;
	dn[6] = (-8.519692053642082*B*C*D*F*G*dbdr) / b2 + (8.519692053642082*B*C*D*F) / b + (8.519692053642082*B*C*D*G) / b + (8.519692053642082*B*C*F*G) / b + (8.519692053642082*B*D*F*G) / b + (8.519692053642082*C*D*F*G) / b;
	dn[8] = (38.28714854897507*A*C*D*F*G*dbdr) / b2 - (38.28714854897507*A*C*D*F) / b - (38.28714854897507*A*C*D*G) / b - (38.28714854897507*A*C*F*G) / b - (38.28714854897507*A*D*F*G) / b - (38.28714854897507*C*D*F*G) / b;
	dn[10] = (-26.495528883381816*A*C*F*E*G*dbdr) / b2 + (26.495528883381816*A*C*F*E) / b + (26.495528883381816*A*C*F*G) / b + (26.495528883381816*A*C*E*G) / b + (26.495528883381816*A*F*E*G) / b + (26.495528883381816*C*F*E*G) / b;
	dn[12] = (-6.543855223902344*A*B*C*E*G*dbdr) / b2 + (6.543855223902344*A*B*C*E) / b + (6.543855223902344*A*B*C*G) / b + (6.543855223902344*A*B*E*G) / b + (6.543855223902344*A*C*E*G) / b + (6.543855223902344*B*C*E*G) / b;
	// Derivadas em relacao a s
	dn[1] = (-26.495528883381816*A*B* C*D*E* dbds) / b2 + (33.22434234515739*A*B* C*D) / b - (12.759574222118465*A*B* C*E) / b - (33.22434234515739*A*B*D* E) / b - (116.08449684783301*A*C*D* E) / b;
	dn[3] = (38.28714854897507*A*B*D* F*E* dbds) / b2 - (48.010565722614885*A*B*D* F) / b - (18.438118967699896*A*B*D* E) / b + (18.438118967699896*A*B*F* E) / b + (167.7469581606903*A*D*F* E) / b;
	dn[5] = (-8.519692053642082*B*D*F* E*G* dbds) / b2 + (37.32720979824761*B*D*F* E) / b + (10.68335592436728*B*D*F* G) / b + (4.102867453090217*B*D* E*G) / b - (4.102867453090217*B*F* E*G) / b - (37.32720979824761*D*F*E* G) / b;
	dn[7] = (-8.519692053642082*B*C* D*F*G* dbds) / b2 + (37.32720979824761*B*C* D*F) / b + (4.102867453090217*B*C* D*G) / b - (4.102867453090217*B*C* F*G) / b - (10.68335592436728*B*D*F* G) / b - (37.32720979824761*C*D*F* G) / b;
	dn[9] = (38.28714854897507*A*C*D* F*G* dbds) / b2 - (167.7469581606903*A*C*D* F) / b - (18.438118967699896*A*C*D* G) / b + (18.438118967699896*A*C*F* G) / b + (48.010565722614885*A*D*F* G) / b;
	dn[11] = (-26.495528883381816*A*C*F* E*G* dbds) / b2 + (116.08449684783301*A*C*F* E) / b + (33.22434234515739*A*C*F* G) / b + (12.759574222118465*A*C* E*G) / b - (33.22434234515739*A*F*E* G) / b;
	dn[13] = (-6.543855223902344*A*B* C*E* G*dbds) / b2 + (28.670503029219354*A*B* C*E) / b + (8.205734906180433*A*B* C*G) / b - (8.205734906180433*A*B* E*G) / b - (28.670503029219354*A*C* E*G) / b;
}