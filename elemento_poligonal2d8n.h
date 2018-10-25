#ifndef ELEMENTO_POLIGONAL2D8N_H
#define ELEMENTO_POLIGONAL2D8N_H

#include "elemento_poligonal2d.h"

/**
@author Fernando Cesar Meira Menandro
*/
class elpol2D8N : public elpol2d
{
private:
	const static int nno = 8;  //Numero de nos
	const static int ptg = 7;
	const static int ptg_tot = ptg;
public:
#ifdef ALEATORIO
	class aleatorio *yg;
#else
	double *yg;
#endif
	elpol2D8N();
	~elpol2D8N();

	int qptg();
	int qnno();
	int qptg_tot();
	void funcao_Forma(double, double, double*, double*);
	void pontos_de_gauss(int, double*, double*, double*);
};

#endif
