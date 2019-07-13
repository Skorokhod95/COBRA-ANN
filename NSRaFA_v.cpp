#include "Algorithms.h"
#include "Functions.h"

#define _USE_MATH_DEFINES

#include<math.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>

using namespace std;


//Норма (евклидова)
double NSRaFA_v::R(double *x, double*y)
{
	int d;
	double sum = 0;
	for (d = 0; d < D; d++)
		sum += pow((x[d] - y[d]), 2.);
	return sqrt(sum);
}

void NSRaFA_v::sort()
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

//Функция инициализации переменных
void NSRaFA_v::initialization()
{
	
	out.open("FA.txt");
	
	int i, d;

	for (int j = 0; j<2; j++)
		G[j] = new double[D];

	nom_x = new double[N];
	nom_pBest = new double[N];
	x = new double*[N];
	for (i = 0; i<N; i++)
		x[i] = new double[D];

	pBest = new double*[N];
	for (i = 0; i<N; i++)
		pBest[i] = new double[D];

	s = new double[D];
	x1 = new double[D];
	x2 = new double[D];
	x3 = new double[D];

	num = new int[N];

	gBest = new double[D];
}

//Функция инициализации популяци объекта
void NSRaFA_v::creation(ANN& ann)
{
	int i, d;

	//Случайная инициаализация положений объектов
	for (i = 0; i < N; i++) {
		for (d = 0; d < D; d++)
			x[i][d] = G[0][d] + ((rand() % 10000) / 10000.)*(G[1][d] - G[0][d]);
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}

	//Определение лучших индивидов для каждого индивида
	for (i = 0; i < N; i++) {
		for (d = 0; d < D; d++)
			pBest[i][d] = x[i][d];
		nom_pBest[i] = nom_x[i];
	}
	// Определение лучшего индивида роя
	for (d = 0; d < D; d++)
		gBest[d] = x[0][d];
	nom_gBest = nom_x[0];
	for (i = 1; i < N; i++)
	if (nom_gBest>nom_x[i]){
		for (d = 0; d < D; d++)
			gBest[d] = x[i][d];
		nom_gBest = nom_x[i];
	}

	//Масштаб каждой расчетной переменной
	for (d = 0; d < D; d++)
		s[d] = G[1][d] - G[0][d];

	k = 3;
}

//Изменение лучших индивидов
void NSRaFA_v::best_decision()
{
	int i, d;

	//Определение лучших индивидов для каждого индивида
	for (i = 0; i < N; i++) {
		if (nom_pBest[i]>nom_x[i]) {
			for (d = 0; d < D; d++)
				pBest[i][d] = x[i][d];
			nom_pBest[i] = nom_x[i];
		}
	}

	// Определение лучшего индивида роя
	if (nom_gBest>nom_x[0]) {
		for (d = 0; d < D; d++)
			gBest[d] = x[0][d];
		nom_gBest = nom_x[0];
	}
}

//Модель случайного привлечения светлячка
void NSRaFA_v::random_attraction(ANN& ann)
{
	int d;

	//Обновление параметров α и β
	ksi = (rand() / double(RAND_MAX)) - 0.5;
	alfa = 0.99*alfa;
	beta = (betaMin + (betaMax - betaMin)*exp(-gamma*pow(R(x[i], x[randF]), 2.)))*((double)t / T);

	//Передвигаем xi к xj, где j=randF
	for (d = 0; d < D; d++)
		x[i][d] = x[i][d] + beta*(x[randF][d] - x[i][d]) + alfa * s[d] * ksi;
	nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
}

// поиск по окресности
void NSRaFA_v::neighborhood_search(ANN& ann)
{
	int i1, i2, i3, i4, d;
	double trans;

	//выбрать два светлячка из к соседского радиуса хi
	if ((i - k) >= 0 && (i + k)<N)
	do {
		i1 = ((double)rand() * ((i + k) - (i - k)) / RAND_MAX + (i - k));
		i2 = ((double)rand() * ((i + k) - (i - k)) / RAND_MAX + (i - k));
	} while (i == i1 || i == i2 || i1 == i2);
	if ((i - k) < 0 && (i + k)<N)
	do {
		i1 = (rand() / double(RAND_MAX))*(i + k);
		i2 = (rand() / double(RAND_MAX))*(i + k);
	} while (i == i1 || i == i2 || i1 == i2);
	if ((i - k) >= 0 && (i + k) >= N)
	do {
		i1 = ((double)rand() * ((N - 1) - (i - k)) / RAND_MAX + (i - k));
		i2 = ((double)rand() * ((N - 1) - (i - k)) / RAND_MAX + (i - k));
	} while (i == i1 || i == i2 || i1 == i2);
	if ((i - k) < 0 && (i + k) >= N)
	do {
		i1 = (rand() / double(RAND_MAX))*(N-1);
		i2 = (rand() / double(RAND_MAX))*(N-1);
	} while (i == i1 || i == i2 || i1 == i2);

	//выбрать два случайных светлячка из всего роя
	do {
		i3 = (rand() / double(RAND_MAX))*(N-1);
		i4 = (rand() / double(RAND_MAX))*(N-1);
	} while (i == i3 || i == i4 || i3 == i4);

	r[0] = (rand() / double(RAND_MAX));
	r[1] = (rand() / double(RAND_MAX))*(1 - r[0]);
	r[2] = 1 - r[0] - r[1];

	r[3] = (rand() / double(RAND_MAX));
	r[4] = (rand() / double(RAND_MAX))*(1 - r[3]);
	r[5] = 1 - r[3] - r[4];

	//Сгенерируем три пробных решения
	for (d = 0; d < D; d++)
	{
		x1[d] = r[0] * x[i][d] + r[1] * pBest[i][d] + r[2] * (x[i1][d] - x[i2][d]); //локальный соседский поисковый оператор
		x2[d] = r[3] * x[i][d] + r[4] * gBest[d] + r[5] * (x[i3][d] - x[i4][d]); // глобальный соседский поисковый оператор
		x3[d] = x[i][d] + Cauchy(); //мутация
	}
	nom_x1 = f(GLOBAL_COUNTER, x1, ann);
	nom_x2 = f(GLOBAL_COUNTER, x2, ann);
	nom_x3 = f(GLOBAL_COUNTER, x3, ann);

	//Выбрать лучшее решение среди х1, х2, х3 и хi
	trans = min(min(nom_x[i], nom_x1), min(nom_x2, nom_x3));
	if (trans == nom_x1) {
		for (d = 0; d < D; d++)
			x[i][d] = x1[d];
		nom_x[i] = nom_x1;
	}
	if (trans == nom_x2) {
		for (d = 0; d < D; d++)
			x[i][d] = x2[d];
		nom_x[i] = nom_x2;
	}
	if (trans == nom_x3) {
		for (d = 0; d < D; d++)
			x[i][d] = x3[d];
		nom_x[i] = nom_x3;
	}

	//t += 2;
}

void NSRaFA_v::NSRaFAv(ANN& ann)
{
		
	for (i = 0; i < N; i++)
	{
		nom_x[i] = f(GLOBAL_COUNTER, x[i], ann);
	}


	sort();
	for (i = 0; i < N; i++)
	{
		//Выбор случайного светлячка
		do {
			randF = (rand() / double(RAND_MAX))*(N-1);
		} while (randF == i);

		//Модель случайного привлечения светлячка
		if (nom_x[randF] < nom_x[i]) {
			random_attraction(ann);
		}
		//Соседский поиск
		else {
			neighborhood_search(ann);
		}
		sort();
		best_decision();
	}
	/*out << '\n' << "Пколение: " << t << '\n';
	for (i = 0; i < N; i++) {
		for (d = 0; d < D; d++)
			out << x[i][d] << '\t';
		out << nom_x[i] << '\n';
	}*/

}

void NSRaFA_v::deleting()
{
	for (i = 0; i<N; i++)
		delete[]x[i];

	for (i = 0; i<N; i++)
		delete[]pBest[i];

	delete[]s;
	delete[]num;
	delete[]x1;
	delete[]x2;
	delete[]x3;

	delete[]gBest;

	delete[]nom_pBest;
	delete[]nom_x;

	for (i = 0; i<2; i++)
		delete[]G[i];
}