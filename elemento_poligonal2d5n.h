#ifndef ELEMENTO_POLIGONAL2D5N_H
#define ELEMENTO_POLIGONAL2D5N_H

#include "elemento_poligonal2d.h"

/**
@author Fernando Cesar Meira Menandro
*/
class elpol2D5N : public elpol2d
{
private:
	const static int nno = 5;  //Numero de nos
	const static int ptg = 16;
	const static int ptg_tot = ptg;
public:
#ifdef ALEATORIO
	class aleatorio *yg;
	void p_processa(aleatorio*);
#else
	double *yg;
	void p_processa(double*);
#endif
	elpol2D5N();
	~elpol2D5N();

	int qptg();
	int qnno();
	int qptg_tot();
	void funcao_Forma(double, double, double*, double*);
	void pontos_de_gauss(int, double*, double*, double*);
	void monta_rigidez();
	void monta_n();
};

#endif
