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
void iba_v::initialization()
{
	out.open("BA.txt");
	
	for (int j = 0; j<2; j++)
		G[j] = new double[D];


	nom_x = new double[N];
	nom_xx = new double[N];
	x = new double*[N];
	for (i = 0; i < N; i++)
		x[i] = new double[D];

	w = new double*[N];
	for (i = 0; i < N; i++)
		w[i] = new double[D];

	xx = new double*[N];
	for (i = 0; i < N; i++)
		xx[i] = new double[D];

	GBest = new double[D];
	wi0 = new double[D];
	wi8 = new double[D];

}

//Функция инициализации популяци объекта
void iba_v::creation(ANN& ann)
{
	int i, d;
	//Случайная инициаализация положений объектов
	for (i = 0; i < N; i++) {
		for (d = 0; d < D; d++)
			x[i][d] = G[0][d] + ((rand() % 10000) / 10000.)*(G[1][d] - G[0][d]);
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}

	for (d = 0; d<D; d++)
		wi0[d] = (G[1][d] - G[0][d]) / 4.;
	for (d = 0; d<D; d++)
		wi8[d] = wi0[d] / 100.;
	a0 = 0.9;
	a8 = 0.6;
	r0 = 0.1;
	r8 = 0.7;
	fMax = 1.0;
	fMin = 0.1;
	sum = 0;

	//Случайная инициаализация значения частот 
	for (i = 0; i<N; i++)
	for (d = 0; d<D; d++)
		w[i][d] = wi0[d] + ((rand() % 10000) / 10000.)*(wi8[d] - wi0[d]);

	//Случайная инициаализация значения громкостей 
	A = a0 + ((rand() % 10000) / 10000.)*(a8 - a0);

	//Случайная инициаализация значения частот повторения импульсов 
	r = r0 + ((rand() % 10000) / 10000.)*(r8 - r0);;


	for (d = 0; d<D; d++)
		GBest[d] = x[0][d];
	nom_GBest = nom_x[0];
	//Лучшее положение поколения
	for (i = 0; i<N; i++) {
		if (nom_x[i] < nom_GBest) {
			for (d = 0; d < D; d++)
				GBest[d] = x[i][d];
			nom_GBest = nom_x[i];
		}
	}

}

void iba_v::sort()
{
	int i, j, k;

	//Сортировка роя 
	for (j = 0; j < N; j++) {
		k = 0;
		for (i = 0; i < N - 1; i++)
		if (nom_x[i]>nom_x[i + 1]) {
			swap(x[i], x[i + 1]);
			swap(nom_x[i], nom_x[i + 1]);
			k++;
		}
		if (k == 0) break;
	}
}

void iba_v::change_of_state(ANN& ann)
{
	int d;
	//Изменение Частоты 
	f1 = fMin + (fMax - fMin)*(rand() / double(RAND_MAX));
	f2 = fMin + (fMax - fMin)*(rand() / double(RAND_MAX));

	//Изменение состояния
	if (f(GLOBAL_COUNTER, x[k], ann) < nom_x[i])
	for (d = 0; d < D; d++)
		x[i][d] = x[i][d] + (GBest[d] - x[i][d])*f1 + (x[k][d] - x[i][d])*f2;
	else
	for (d = 0; d < D; d++)
		x[i][d] = x[i][d] + (GBest[d] - x[i][d])*f1;
	nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
}

void iba_v::local_search(ANN& ann)
{
	int d;
	double q;
	q = rand() / double(RAND_MAX);
	if (q>r)
	for (d = 0; d < D; d++)
	{
		xx[i][d] = x[i][d] + aMid*(((double)rand() * (1 - (-1)) / RAND_MAX + (-1)))*w[i][d];
		w[i][d] = ((wi0[d] - wi8[d]) / (1 - T))*(t - T) + wi8[d];
	}
	nom_xx[i] = f(GLOBAL_COUNTER, xx[i], ann);
}

void iba_v::best_decision_global()
{
	int ii, dd;
	//Лучшее положение поколения
	if (nom_x[0] < nom_GBest) {
		for (dd = 0; dd < D; dd++)
			GBest[dd] = x[0][dd];
		nom_GBest = nom_x[0];
	}
}

void iba_v::accept_local_search()
{
	int d;
	double q;
	q = rand() / double(RAND_MAX);
	if (q < A && nom_GBest < nom_x[i])
	{
		for (d = 0; d < D; d++)
			x[i][d] = xx[i][d];
		nom_x[i] = nom_xx[i];
		r = ((r0 - r8) / (1 - T))*(t - T) + r8;
		A = ((a0 - a8) / (1 - T))*(t - T) + a8;
	}
}

void iba_v::Ibav(ANN& ann)
{
	for (i = 0; i < N; i++)
	{
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}

	sort();
	sum += A;
	aMid = (double)A / (t + 1);

	for (i = 0; i < N; i++)
	{
		k = -1;
		//случайная мышь
		do {
			k1 = rand() % N;
			if (k1 != i) k = k1;
		} while (k == -1);

		change_of_state(ann);
		local_search(ann);
		this->accept_local_search();
	}
	sort();
	best_decision_global();

	/*out << '\n' << "Пколение: " << t << '\n';
	for (i = 0; i < N; i++) {
		for (d = 0; d < D; d++)
			out << x[i][d] << '\t';
		out << nom_x[i] << '\n';
	}*/
}

void iba_v::deleting()
{
	for (i = 0; i<N; i++)
		delete[]x[i];

	for (i = 0; i<N; i++)
		delete[]xx[i];

	for (i = 0; i<N; i++)
		delete[]w[i];

	delete[]GBest;


	delete[]nom_x;
	delete[]nom_xx;

	for (i = 0; i<2; i++)
		delete[]G[i];
}
