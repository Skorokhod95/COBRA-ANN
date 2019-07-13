#include "Algorithms.h"
#include "Functions.h"

#define _USE_MATH_DEFINES

#include<math.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>

using namespace std;



double CSKH_v::K_norm(int i, double X)
{
	return ((nom_x[i] - X) / (nom_gWorst - nom_gBest));
}

double CSKH_v::X_norm(int i, double *X, int D)
{
	int d;
	double sum = 0, normXij;

	for (d = 0; d < DP; d++)
		sum += pow((x[i][d] - X[d]), 2.);
	normXij = sqrt(sum);
	return (x[i][D] - X[D]) / (normXij + 0.0001);
}

//Функция инициализации переменных
void CSKH_v::initialization()
{
	out.open("CSA.txt");
	int i, d;

	nom_x = new double[NP];
	nom_x_Neighbors = new double[NP];
	nom_phantom1 = new double[NP];
	nom_P1 = new double[NP];
	nom_P2 = new double[NP];

	phantom = new double[DP];
	phantom1 = new double*[NP];
	for (i = 0; i<NP; i++)
		phantom1[i] = new double[DP];

	K = new bool*[NP];
	for (i = 0; i<NP; i++)
		K[i] = new bool[DP];

	x = new double*[NP];
	for (i = 0; i<NP; i++)
		x[i] = new double[DP];

	xNeighbors = new double*[NP];
	for (i = 0; i<NP; i++)
		xNeighbors[i] = new double[DP];
	norm = new double[NP];
	ds = new double[NP];

	N = new double[DP];
	Nold = new double*[NP];
	for (i = 0; i<NP; i++)
		Nold[i] = new double[DP];

	F = new double[DP];
	Fold = new double*[NP];
	for (i = 0; i<NP; i++)
		Fold[i] = new double[DP];

	D = new double[DP];

	beta = new double[DP];
	betaBest = new double[DP];
	betaFood = new double[DP];

	alfa = new double[DP];
	alfaLocal = new double[DP];
	alfaTarget = new double[DP];

	xfood = new double[DP];

	gBest = new double[DP];
	gWorst = new double[DP];
	num = new int[NP];

	step = new double[DP];
	P1 = new double*[NP];
	for (i = 0; i<NP; i++)
		P1[i] = new double[DP];
	P2 = new double*[NP];
	for (i = 0; i<NP; i++)
		P2[i] = new double[DP];

	for (int j = 0; j<2; j++)
		G[j] = new double[DP];
}

//Функция инициализации популяци объекта
void CSKH_v::creation(ANN& ann)
{
	int i, d;
	double sum = 0;

	//Случайная инициаализация положений объектов
	for (i = 0; i < NP; i++) {
		for (d = 0; d < DP; d++)
			x[i][d] = G[0][d] + ((rand() % 10000) / 10000.)*(G[1][d] - G[0][d]);
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}

	// Определение лучшего индивида роя
	for (d = 0; d < DP; d++)
		gBest[d] = x[0][d];
	nom_gBest = f(GLOBAL_COUNTER, gBest, ann);
	for (i = 1; i < NP; i++)
		if (nom_gBest>nom_x[i]) {
			for (d = 0; d < DP; d++)
				gBest[d] = x[i][d];
			nom_gBest = nom_x[i];
		}
	// Определение худшего индивида роя
	for (d = 0; d < DP; d++)
		gWorst[d] = x[0][d];
	nom_gWorst = f(GLOBAL_COUNTER, gWorst, ann);
	for (i = 1; i < NP; i++)
		if (nom_gWorst < nom_x[i]) {
			for (d = 0; d < DP; d++)
				gWorst[d] = x[i][d];
			nom_gWorst = nom_x[i];
		}

	Vf = 0.02;
	DMax = 0.005;
	NMax = 0.01;
	pa = 0.1;
	omegaF = 0.9;
	omegaN = omegaF;

	Ct = 0.5;	//0..2
	for (d = 0; d < DP; d++)
		sum += (G[1][d] - G[0][d]);
	dt = Ct*sum;

	lyamda = 2; //1..3
	B = 0.2;

	for (i = 0; i < NP; i++)
	for (d = 0; d < DP; d++)
	{
		Nold[i][d] = 0;
		Fold[i][d] = 0;
	}
}

//Сортировка роя
void CSKH_v::sort()
{
	int i, j, k;

	//Сортировка роя 
	for (j = 0; j < NP; j++) {
		k = 0;
		for (i = 0; i < NP - 1; i++)
		if (nom_x[i]>nom_x[i+1]) {
			swap(x[i], x[i + 1]);
			swap(nom_x[i], nom_x[i + 1]);
			k++;
		}
		if (k == 0) break;
	}
}

//Запись и перестановка лучших индивидов
void CSKH_v::work_with_best(char choise)
{
	int i, d;

	if (choise == 's')
	{
		for (i = 0; i < KEEP; i++)
		for (d = 0; d < DP; d++)
			phantom1[i][d] = x[i][d];
	}
	else
	{
		for (i = 0; i < KEEP; i++)
		for (d = 0; d < DP; d++)
			x[NP - 1 - i][d] = phantom1[i][d];
	}
}

//Алгоритм KU
void CSKH_v::CU(ANN& ann)
{
	int d, j;

	//Поиск состояния (решения) через алгоритм кукушки (случайных полетов Леви)
	for (d = 0; d < DP; d++)
		phantom[d] = x[i][d] + B*Levy(lyamda);

	j = (rand() / double(RAND_MAX))*(NP - 1);
	if (nom_x[j]>nom_x[i]) {
		for (d = 0; d < DP; d++)
			x[j][d] = phantom[d];
		nom_x[j] = nom_x[i];
		}
	else
		KH(ann);
}

//Перемена состояния через алгоритм криля
void CSKH_v::KH(ANN& ann)
{
	int d, j;
	double sum1 = 0, sum2 = 0;
	double nom_Xfood;

	Cbest = 2 * ((rand() / double(RAND_MAX)) + ((double)t / T));
	Cfood = 2 * (1 - ((double)t / T));

	for (d = 0; d < DP; d++) {
		for (j = 0; j < NP; j++)
		{
			sum1 += (x[j][d] / nom_x[j]);
			sum2 += (1 / nom_x[j]);
		}
		xfood[d] = sum1 / sum2;
		sum1 = 0;
		sum2 = 0;
	}
	nom_Xfood = f(GLOBAL_COUNTER, xfood, ann);

	neighborhood_search();

	for (d = 0; d < DP; d++)
	{
		betaBest[d] = K_norm(i, nom_gBest)*X_norm(i,gBest, d);
		betaFood[d] = Cfood*K_norm(i, nom_Xfood)*X_norm(i, xfood, d);
		alfaTarget[d] = Cbest*K_norm(i, nom_gBest)*X_norm(i, gBest, d);
		for (j = 0; j < NN; j++)
			sum1 += (K_norm(i, nom_x_Neighbors[j])*X_norm(i, xNeighbors[j], d));
		alfaLocal[d] = sum1;
		sum1 = 0;

		alfa[d] = alfaLocal[d] + alfaTarget[d];
		beta[d] = betaBest[d] + betaFood[d];

		N[d] = NMax*alfa[d] + omegaN*Nold[i][d];
		F[d] = Vf*beta[d] + omegaF*Fold[i][d];
		D[d] = DMax*((double)rand() * (1 - (-1)) / RAND_MAX + (-1));

		x[i][d] = x[i][d] + dt*(N[d] + F[d] + D[d]);

		Fold[i][d] = F[d];
		Nold[i][d] = N[d];
	}
	nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
}

//Изменение лучших индивидов
void CSKH_v::best_decision()
{
	int i, d;

	//Определение худших индивидов 
	for (i = 1; i < NP; i++)
	if (nom_gWorst < nom_x[i]){
		for (d = 0; d < DP; d++)
			gWorst[d] = x[i][d];
		nom_gWorst = nom_x[i];
	}

	// Определение лучшего индивида роя
	for (i = 0; i < NP; i++)
	if (nom_gBest>nom_x[i]) {
		for (d = 0; d < DP; d++)
			gBest[d] = x[i][d];
		nom_gBest = nom_x[i];
	}
}

//Поиск по окрестностям
void CSKH_v::neighborhood_search()
{
	int j, jj = 0, d;
	double sum = 0, sum2 = 0;

	for (j = 0; j < NP; j++) {
		for (d = 0; d < DP; d++)
			sum2 += pow((x[i][d] - x[j][d]), 2.);
		norm[j] = sqrt(sum2);
		sum += norm[j];
		sum2 = 0;
	}
	ds[i] = sum / ((double)2 * NP);
	sum = 0;

	for (j = 0; j < NP; j++)
	if (norm[j] <= ds[i] && i != j)
	{
		for (d = 0; d < DP; d++)
			xNeighbors[jj][d] = x[j][d];
		nom_x_Neighbors[jj] = nom_x[i];
		jj++;
	}
	NN = jj;
	jj = 0;
}

//Алгоритм CA
void CSKH_v::CA(ANN& ann)
{
	int i, d, random;

	for (i = 0; i < NP; i++) {
		for (d = 0; d < DP; d++)
		{
			P1[i][d] = x[i][d];
			P2[i][d] = x[i][d];
		}
		nom_P1[i] = nom_x[i];
		nom_P2[i] = nom_x[i];
	}

	for (i = NP - 1; i >= 0; i--) {
		random = (rand() / double(RAND_MAX))*i;
		swap(P1[i], P1[random]);
		swap(nom_P1[i], nom_P1[random]);
	}
	for (i = NP - 1; i >= 0; i--) {
		random = (rand() / double(RAND_MAX))*i;
		swap(P2[i], P2[random]);
		swap(nom_P2[i], nom_P2[random]);
	}

	for (i = 0; i < NP; i++)
	{
		for (d = 0; d < DP; d++)
		{
			if ((rand() / double(RAND_MAX))>pa)
				K[i][d] = true;
			else
				K[i][d] = false;
		}
	}

	for (i = 0; i < NP; i++)
	{
		for (d = 0; d < DP; d++)
		{
			step[d] = (rand() / double(RAND_MAX)) * (P1[i][d] - P2[i][d]);
			phantom[d] = x[i][d] + step[d] * K[i][d];
		}
		nom_phantom = f(GLOBAL_COUNTER, phantom, ann);
		if (nom_phantom < nom_x[i]) {
			for (d = 0; d < DP; d++)
				x[i][d] = phantom[d];
			nom_x[i] = nom_phantom;
		}
	}

}

void CSKH_v::CSKHv(ANN& ann)
{
	for (i = 0; i < NP; i++)
	{
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}
		
		sort();
		work_with_best('s');
		for (i = 0; i < NP; i++)
			CU(ann);
		CA(ann);
		work_with_best('r');
		sort();
		best_decision();

		/*out << '\n' << "Пколение: " << t << '\n';
		for (i = 0; i < NP; i++) {
			for (d = 0; d < DP; d++)
				out << x[i][d] << '\t';
			out << nom_x[i] << '\n';
		}*/
}

void CSKH_v::deleting()
{
	delete[]nom_x;
	delete[]nom_x_Neighbors;
	delete[]nom_phantom1;
	delete[]nom_P1;
	delete[]nom_P2;
	
	delete[]xfood;

	delete[]num;

	for (i = 0; i<NP; i++)
		delete[]K[i];

	for (i = 0; i<NP; i++)
		delete[]x[i];

	for (i = 0; i<NP; i++)
		delete[]xNeighbors[i];
	delete[]norm;
	delete[]ds;


	delete[]F;
	for (i = 0; i<NP; i++)
		delete[]Fold[i];

	delete[]D;

	delete[]phantom;
	for (i = 0; i<NP; i++)
		delete[]phantom1[i];

	delete[]beta;
	delete[]betaBest;
	delete[]betaFood;

	delete[]alfa;
	delete[]alfaLocal;
	delete[]alfaTarget;

	delete[]N;
	for (i = 0; i<NP; i++)
		delete[]Nold[i];

	delete[]gBest;

	delete[]gWorst;

	delete[]step;
	for (i = 0; i<NP; i++)
		delete[]P1[i];
	for (i = 0; i<NP; i++)
		delete[]P2[i];

	for (i = 0; i<2; i++)
		delete[]G[i];
}