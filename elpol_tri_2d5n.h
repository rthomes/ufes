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

#ifndef ELPOL_TRI_2D5N_H
#define ELPOL_TRI_2D5N_H

#include "elpol_tri_2d.h"

class elpol_tri_2D5N : public elpol_tri_2d
{
private:
	const static int nno = 5;  //Numero de nos
	//const static int ptg = 3;  //Numero de Pontos de Gauss nas diferentes direcoes
public:
    #ifdef ALEATORIO
   class aleatorio *yg;
#else
   double *yg;
#endif
	elpol_tri_2D5N();
	~elpol_tri_2D5N();

	int qnno();
	void funcao_Forma(double, double, double*, double*);
};

#endif
