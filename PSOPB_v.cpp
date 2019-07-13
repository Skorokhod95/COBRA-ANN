#include "Algorithms.h"
#include "Functions.h"

#define _USE_MATH_DEFINES

#include<math.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>

using namespace std;

//Функция инициализации переменных
void PSOPB_v::initialization()
{
	out.open("PSO.txt");

	Np = N / 2;
	Nh = N - Np;

	swarmP = new double*[Np];
	for (i = 0; i<Np; i++)
		swarmP[i] = new double[D];
	nom_swarmP = new double[Np];


	swarmH = new double*[Nh];
	for (i = 0; i<Nh; i++)
		swarmH[i] = new double[D];
	nom_swarmH = new double[Nh];


	VH = new double*[Nh];
	for (i = 0; i<Nh; i++)
		VH[i] = new double[D];

	VP = new double*[Np];
	for (i = 0; i<Np; i++)
		VP[i] = new double[D];

	piH = new double*[Nh];
	for (i = 0; i<Nh; i++)
		piH[i] = new double[D];

	piP = new double*[Np];
	for (i = 0; i<Np; i++)
		piP[i] = new double[D];

	pgP = new double*[T];
	for (i = 0; i<T; i++)
		pgP[i] = new double[D];

	pgH = new double[D];
	vMax = new double[D];

	numP = new int[Np];
	numH = new int[Nh];

	pc = new double[D];
	pn = new double[D];

	GB = new double[D];


	nom_piH = new double[Nh];

	nom_piP = new double[Np];

	nom_pgP = new double[T];

	for (int j = 0; j<2; j++)
		G[j] = new double[D];

}

//Функция инициализации популяци объекта
void PSOPB_v::creation(ANN& ann)
{
	Np = N / 2;
	Nh = N - Np;
	lyambda = 0.729;
	c1 = 2.05; c2 = c1;
	c11 = 1.367; c12 = c11; c13 = c11;
	w = 0.7968;
	gamma = 0.5;
	ita = 7;
	pcMax = 0.02; pcMin = 0.005;
	ksiMax = 0.5; ksiMin = 0.3;
	dMax = 1; dMin = 0.1;

	//Случайная инициаализация положений объектов
	for (i = 0; i < Nh; i++) {
		for (d = 0; d < D; d++)
			swarmH[i][d] = G[0][d] + ((rand() % 10000) / 10000.)*(G[1][d] - G[0][d]);
		nom_swarmH[i] = f(GLOBAL_COUNTER, swarmH[i], ann);
		numH[i] = i;
	}
	for (i = 0; i < Np; i++) {
		for (d = 0; d < D; d++)
			swarmP[i][d] = G[0][d] + ((rand() % 10000) / 10000.)*(G[1][d] - G[0][d]);
		nom_swarmP[i] = f(GLOBAL_COUNTER, swarmP[i], ann);
		numP[i] = i;
	}

	//Инициализация начальных скоростей объектов
	for (i = 0; i<Nh; i++)
	for (d = 0; d<D; d++)
		VH[i][d] = 0;
	for (i = 0; i<Np; i++)
	for (d = 0; d<D; d++)
		VP[i][d] = 0;

	for (d = 0; d < D; d++)
		vMax[d] = 0.25 * (G[1][d] - G[0][d]);

	//Инициализация лучших положений i-ой точки
	for (i = 0; i < Nh; i++) {
		for (d = 0; d < D; d++)
			piH[i][d] = swarmH[i][d];
		nom_piH[i] = nom_swarmH[i];
	}
	for (i = 0; i < Np; i++) {
		for (d = 0; d<D; d++)
			piP[i][d] = swarmP[i][d];
		nom_piP[i] = nom_swarmP[i];
	}

	for (d = 0; d<D; d++)
		pgH[d] = swarmH[0][d];
	nom_pgH = nom_swarmH[0];
	//Лучшее положение роя хозяев
	for (i = 0; i<Nh; i++) {
		if (nom_swarmH[i] < nom_pgH) {
			for (d = 0; d < D; d++)
				pgH[d] = swarmH[i][d];
			nom_pgH = nom_swarmH[i];
		}
	}


	for (d = 0; d<D; d++)
		pgP[0][d] = swarmP[0][d];
	nom_pgP[0] = nom_swarmP[0];
	//Лучшее положение роя паразитов
	for (i = 1; i<Np; i++) {
		if (nom_swarmP[i] < nom_pgP[0]) {
			for (d = 0; d < D; d++)
				pgP[0][d] = swarmP[i][d];
			nom_pgP[0] = nom_swarmP[i];
		}
	}


	//Лучшее историческое положение
	if (nom_pgH < nom_pgP[0]) {
		for (d = 0; d < D; d++)
			GB[d] = pgH[d];
		nom_GB = nom_pgH;
	}
	else {
		for (d = 0; d < D; d++)
			GB[d] = pgP[0][d];
		nom_GB = nom_pgP[0];
	}

}

//Функция изменения состояния объекта: скорости и координаты
void PSOPB_v::change_of_state(ANN& ann)
{
	//Обновление состояния роя паразитов
	for (i = 0; i < Np; i++) {
		for (d = 0; d < D; d++)
		{
			VP[i][d] = lyambda * (w*VP[i][d] + c1*(rand() / double(RAND_MAX))*(piP[i][d] - swarmP[i][d]) + c2*(rand() / double(RAND_MAX))*(pgP[t][d] - swarmP[i][d]));
			swarmP[i][d] = VP[i][d] + swarmP[i][d];
		}
		nom_swarmP[i] = f(GLOBAL_COUNTER, swarmP[i], ann);
	}

	//Обновление состояния роя паразитов
	for (i = 0; i < Nh; i++) {
		if (nom_pgH >= nom_pgP[t])
		for (d = 0; d < D; d++)
			VH[i][d] = lyambda * (VH[i][d] + c11*(rand() / double(RAND_MAX))*(piH[i][d] - swarmH[i][d]) + c12*(rand() / double(RAND_MAX))*(pgH[d] - swarmH[i][d]) + c13*(rand() / double(RAND_MAX))*(pgP[t][d] - swarmH[i][d]));
		else
		for (d = 0; d < D; d++)
			VH[i][d] = lyambda * (VH[i][d] + c1*(rand() / double(RAND_MAX))*(piH[i][d] - swarmH[i][d]) + c2*(rand() / double(RAND_MAX))*(pgH[d] - swarmH[i][d]));
		for (d = 0; d < D; d++)
			swarmH[i][d] = VH[i][d] + swarmH[i][d];
		nom_swarmH[i] = f(GLOBAL_COUNTER, swarmH[i], ann);
	}
}

//Функция поиска лучшего решения для каждой точки, для всего поколения
void PSOPB_v::best_decisions()
{
	//Инициализация лучших положений i-ой точки
	for (i = 0; i < Nh; i++)
	if (nom_swarmH[i] < nom_piH[i]) {
		for (d = 0; d < D; d++)
			piH[i][d] = swarmH[i][d];
		nom_piH[i] = nom_swarmH[i];
	}
	for (i = 0; i < Np; i++)
	if (nom_swarmP[i] < nom_piP[i]) {
		for (d = 0; d < D; d++)
			piP[i][d] = swarmP[i][d];
		nom_piP[i] = nom_swarmP[i];
	}

	//Лучшее положение роя хозяев
	if (nom_pgH>nom_swarmH[0]) {
		for (d = 0; d < D; d++)
			pgH[d] = swarmH[0][d];
		nom_pgH = nom_swarmH[0];
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Лучшее положение роя паразитов
	if (t >= 1) {
		for (d = 0; d < D; d++)
			pgP[t][d] = pgP[t - 1][d];
		nom_pgP[t] = nom_pgP[t - 1];
		if (nom_swarmP[0] < nom_pgP[t]) {
			for (d = 0; d < D; d++)
				pgP[t][d] = swarmP[0][d];
			nom_pgP[t] = nom_swarmP[0];
		}
	}
}

//Лучшее историческое положение
void PSOPB_v::best_decision()
{
	if (nom_pgP[t]<nom_GB || nom_pgH<nom_GB) {

		if (nom_pgH < nom_pgP[t]) {
			for (d = 0; d < D; d++)
				GB[d] = pgH[d];
			nom_GB = nom_pgH;
		}
		else {
			for (d = 0; d < D; d++)
				GB[d] = pgP[t][d];
			nom_GB = nom_pgP[t];
		}
	}
}

//Сортировка обоих роев
void PSOPB_v::sort()
{
	int i, j, k, num;

	//Сортировка роя хозяев
	for (j = 0; j < Nh; j++) {
		k = 0;
		for (i = 0; i < Nh - 1; i++)
		if (nom_swarmH[i]>nom_swarmH[i + 1]) {
			swap(swarmH[i], swarmH[i + 1]);
			swap(nom_swarmH[i], nom_swarmH[i + 1]);
			k++;
		}
		if (k == 0) break;
	}

	//Сортировка роя паразитов
	for (j = 0; j < Np; j++) {
		k = 0;
		for (i = 0; i < Np - 1; i++)
		if (nom_swarmP[i]>nom_swarmP[i + 1]) {
			swap(swarmP[i], swarmP[i + 1]);
			swap(nom_swarmP[i], nom_swarmP[i + 1]);
			k++;
		}
		if (k == 0) break;
	}
}

//Обмен лучшими и худшими индивидами между роем хозяев и роем паразитов
void PSOPB_v::exchange()
{
	double **phantomH, **phantomP;
	double *nom_phantomH, *nom_phantomP;

	nom_phantomP = new double[Np];
	nom_phantomH = new double[Nh];
	phantomP = new double*[Np];
	for (i = 0; i<Np; i++)
		phantomP[i] = new double[D];
	phantomH = new double*[Nh];
	for (i = 0; i<Nh; i++)
		phantomH[i] = new double[D];

	for (i = 0; i < Np; i++) {
		for (d = 0; d < D; d++)
			phantomP[i][d] = swarmP[i][d];
		nom_phantomP[i] = nom_swarmP[i];
	}
	for (i = 0; i < Nh; i++) {
		for (d = 0; d < D; d++)
			phantomH[i][d] = swarmH[i][d];
		nom_phantomH[i] = nom_swarmH[i];
	}

	//Рой паразитов получает ξ*Nh лучших индивидов роя хозяев
	for (i = 0; i < Nh*ksi; i++) {
		for (d = 0; d < D; d++)
			swarmP[i][d] = phantomH[i][d];
		nom_swarmP[i] = nom_phantomH[i];
	}
	//Рой хозяев получает ξ*Nh худших индивидов роя паразитов
	for (i = 0; i < Nh*ksi; i++) {
		for (d = 0; d < D; d++)
			swarmH[i][d] = phantomP[Np - 1 - i][d];
		nom_swarmH[i] = nom_phantomP[i];
	}

	delete[]nom_phantomH;
	delete[]nom_phantomP;
	for (i = 0; i<Nh; i++)
		delete[]phantomH[i];
	for (i = 0; i<Np; i++)
		delete[]phantomP[i];
}

//Удаление худших индивидов из роев хозяев
void PSOPB_v::reinitialization(ANN& ann)
{
	psi = 0.3 - 0.25*((double)t / T);

	for (i = 0; i < gamma*Nh; i++)
	{
		if ((rand() / double(RAND_MAX)) <= 0.5) {
			for (d = 0; d < D; d++)
				swarmH[Nh - 1 - i][d] = G[0][d] + (G[1][d] - G[0][d])*(rand() / double(RAND_MAX));
			nom_swarmH[Nh - 1 - i] = f(GLOBAL_COUNTER, swarmH[Nh - 1 - i], ann);
		}
		else {
			for (d = 0; d < D; d++)
				swarmH[Nh - 1 - i][d] = GB[d] + psi * (G[1][d] - G[0][d]) * (2 * (rand() / double(RAND_MAX)) - 1);
			nom_swarmH[Nh - 1 - i] = f(GLOBAL_COUNTER, swarmH[Nh - 1 - i], ann);
		}
	}
}

//Мутация случайного индивида через скорость
void PSOPB_v::mutation(ANN& ann)
{
	//Выбираем случайного индивида из роя паразитов для его мутации 
	int randP;
	randP = (rand() / double(RAND_MAX)) * (Np - 1);

	//"Мутируем" скорость
	for (d = 0; d < D; d++) {
		if ((rand() / double(RAND_MAX)) < 0.5)
			VP[randP][d] = 0.5*vMax[d] * (rand() / double(RAND_MAX));
		else
			VP[randP][d] = -0.5*vMax[d] * (rand() / double(RAND_MAX));
		//Смена положения
		swarmP[randP][d] = swarmP[randP][d] + VP[randP][d];
	}
	nom_swarmH[randP] = f(GLOBAL_COUNTER, swarmH[randP], ann);
}

//Мутация случайного гена
void PSOPB_v::positional_mutation(ANN& ann)
{
	int d;
	double nom_pc, nom_pn;
	d_mutation = (rand() / double(RAND_MAX))*(D - 1);

	for (d = 0; d < D; d++)
	{
		if (d == d_mutation) {
			pn[d_mutation] = pgP[t][d_mutation] + (G[1][d_mutation] - G[0][d_mutation])*Gaussian(0, 1);
			pc[d_mutation] = pgP[t][d_mutation] + Cauchy();

			//Обработка граничных ограничений
			if (pn[d_mutation] < G[0][d_mutation]) pn[d_mutation] = min(G[1][d_mutation], 2 * G[0][d_mutation] - pn[d_mutation]);
			if (pn[d_mutation] > G[1][d_mutation]) pn[d_mutation] = max(G[0][d_mutation], 2 * G[1][d_mutation] - pn[d_mutation]);

			if (pc[d_mutation] < G[0][d_mutation]) pc[d_mutation] = min(G[1][d_mutation], 2 * G[0][d_mutation] - pc[d_mutation]);
			if (pc[d_mutation] > G[1][d_mutation]) pc[d_mutation] = max(G[0][d_mutation], 2 * G[1][d_mutation] - pc[d_mutation]);
		}
		else {
			pc[d] = pgP[t][d];
			pn[d] = pgP[t][d];
		}
	}
	nom_pc = f(GLOBAL_COUNTER, pc, ann);
	nom_pn = f(GLOBAL_COUNTER, pn, ann);

	//смена минимума
	if (nom_pn > nom_pc)
	if (nom_pgP[t] > nom_pc) {
		for (d = 0; d < D; d++)
			pgP[t][d] = pc[d];
		nom_pgP[t] = nom_pc;
	}
	else
	if (nom_pgP[t]>nom_pn) {
		for (d = 0; d < D; d++)
			pgP[t][d] = pn[d];
		nom_pgP[t] = nom_pn;
	}
}

void PSOPB_v::PSOPBv(ANN& ann)
{
	for ( i = 0; i < Nh; i++)
	{
		nom_swarmH[i] = f(GLOBAL_COUNTER, swarmH[i], ann);
	}
	for (i = 0; i < Np; i++)
	{
		nom_swarmP[i] = f(GLOBAL_COUNTER, swarmP[i], ann);
	}
	sort();
	best_decisions();
	change_of_state(ann);
	Pc = pcMin + (pcMax - pcMin) * pow(((double)t / T), 3.);
	if ((rand() / double(RAND_MAX)) <= Pc)
	{
		sort();
		ksi = ksiMax - (ksiMax - ksiMin)*((double)t / T);
		exchange();
		reinitialization(ann);
	}
	if (t >= ita)
	if (nom_pgP[t] - nom_pgP[t - ita] == 0)
		mutation(ann);
	positional_mutation(ann);
	sort();
	best_decisions();
	best_decision();

	/*out << '\n' << "Пколение: " << t << '\t' << N << '\n';
	for (i = 0; i < Nh; i++) {
		for (d = 0; d < D; d++)
			out << swarmH[i][d] << '\t';
		out << nom_swarmH[i] << '\n';
	}
	for (i = 0; i < Np; i++) {
		for (d = 0; d < D; d++)
			out << swarmP[i][d] << '\t';
		out << nom_swarmP[i] << '\n';
	}*/

}

void PSOPB_v::deleting()
{
	for (i = 0; i<Nh; i++)
		delete[]swarmH[i];
	for (i = 0; i<Np; i++)
		delete[]swarmP[i];

	for (i = 0; i<Nh; i++)
		delete[]VH[i];
	for (i = 0; i<Np; i++)
		delete[]VP[i];

	for (i = 0; i<Nh; i++)
		delete[]piH[i];
	for (i = 0; i<Np; i++)
		delete[]piP[i];

	for (i = 0; i<T; i++)
		delete[]pgP[i];

	delete[]pgH;
	delete[]vMax;

	delete[]numH;
	delete[]numP;

	delete[]pc;
	delete[]pn;

	delete[]GB;

	delete[]nom_pgP;
	delete[]nom_piH;
	delete[]nom_piP;
	delete[]nom_swarmH;
	delete[]nom_swarmP;

	for (i = 0; i<2; i++)
		delete[]G[i];
}