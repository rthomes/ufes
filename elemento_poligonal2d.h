#ifndef ELEMENTO_POLIGONAL2D_H
#define ELEMENTO_POLIGONAL2D_H

#include "isop2d.h"

class elpol2d : public isop2d
{
private:
	//const static int ptg = 25;
public:
#ifdef ALEATORIO
	class aleatorio *yg;
	void p_processa(aleatorio*);
#else
	double *yg;
	void p_processa(double*);
#endif
	double *rpg, *spg, *wpg;

	elpol2d();
	elpol2d(int, int, int);
	~elpol2d();

	virtual int qptg() = 0;
	virtual int qptg_tot() = 0; // Adicionado por Renan
	virtual void pontos_de_gauss(int, double*, double*, double*);
	virtual void monta_rigidez();
	virtual void monta_n();
	virtual void funcao_Forma(double r, double s, double *N, double *dn) = 0;
	int tri;
};

#endif
