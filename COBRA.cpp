#include "Algorithms.h"
#include "Functions.h"

#define _USE_MATH_DEFINES

#include<math.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <omp.h>

using namespace std;


void COBRA::input(PSOPB_v & algorythm1, iba_v & algorythm2, NSRaFA_v& algorythm3, CSKH_v& algorythm4, ANN& ann)
{
	int j, d, i;


	//Ввод размерности вектора
	//cout << "Введите размерность задачи: ";
	//cin >> D;
	algorythm1.D = D;
	algorythm2.D = D;
	algorythm3.D = D;
	algorythm4.DP = D;

	algorythm1.N = N_MAX;
	algorythm2.N = N_MAX;
	algorythm3.N = N_MAX;
	algorythm4.NP = N_MAX;

	//cout << '\n' << "Введите максимальное поколение: ";
	//cin >> T;
	T = 100;
	//Количество поколений
	algorythm1.T = T;
	algorythm2.T = T;
	algorythm3.T = T;
	algorythm4.T = T;

	for (i = 0; i < M; i++) {
		Gbest[i] = new double[D];
	}
	GBest = new double[D];


	algorythm1.initialization();
	algorythm2.initialization();
	algorythm3.initialization();
	algorythm4.initialization();

	
	//Ввод границ поиска
	for (j = 0; j<2; j++)
		G[j] = new double[D];
	//cout << '\n' << "Введите границы поиска: ";
	for (d = 0; d < D; d++) {
		G[0][d] = -1;
		G[1][d] = 1;
	}
	for (j = 0; j < 2; j++) {
		for (d = 0; d < D; d++) {
			algorythm1.G[j][d] = G[j][d];
			algorythm2.G[j][d] = G[j][d];
			algorythm3.G[j][d] = G[j][d];
			algorythm4.G[j][d] = G[j][d];
		}
	}

	algorythm1.creation(ann);
	algorythm2.creation(ann);
	algorythm3.creation(ann);
	algorythm4.creation(ann);
	
	for (d = 0; d < D; d++)
		GBest[d] = ((double)rand() * (G[1][d] - G[0][d]) / RAND_MAX + G[0][d]);
	nom_GBest = f(GLOBAL_COUNTER, GBest, ann);

	

	algorythm1.N = N_MIN;
	algorythm2.N = N_MIN;
	algorythm3.N = N_MIN;
	algorythm4.NP = N_MIN;

		
	for (j = 0; j < M; j++) {
		nom_x[j] = new double[N_MAX];
		x[j] = new double*[N_MAX];
		for (i = 0; i < N_MAX; i++)
			x[j][i] = new double[D];
	}

	for (i = 0; i < M; i++)
		N[i] = N_MIN;
}

void COBRA::work_of_algorythms(PSOPB_v & algorythm1, iba_v & algorythm2, NSRaFA_v& algorythm3, CSKH_v& algorythm4, ANN& ann)
{
	int i1,i2,i3,i4, j,d1,d2,d3,d4,t1,t2,t3,t4,i,d;
	//ofstream fout;
//	fout.open("out2.txt");
	t1 = t; t2 = t; t3 = t; t4 = t;

	//Размерность популяции
	algorythm1.N = N[0];
	algorythm1.Nh = N[0] / (double)2;
	algorythm1.Np = N[0] - algorythm1.Nh;
				if (t1 >= 1) {
					for (i1 = 0; i1 < N[0]; i1++)
					{
						for (d1 = 0; d1 < D; d1++) {
							if (i1 < algorythm1.Nh)
								algorythm1.swarmH[i1][d1] = x[0][i1][d1];
							else algorythm1.swarmP[i1 - algorythm1.Nh][d1] = x[0][i1][d1];
						}
					}
				}
				algorythm1.PSOPBv(ann);

				for (i1 = 0; i1 < N[0]; i1++)
				{
					for (d1 = 0; d1 < D; d1++) {
						if (i1 < algorythm1.Nh) x[0][i1][d1] = algorythm1.swarmH[i1][d1];
						else x[0][i1][d1] = algorythm1.swarmP[i1 - algorythm1.Nh][d1];
					}

				}
				for (d1 = 0; d1 < D; d1++)
					Gbest[0][d1] = algorythm1.GB[d1];
				nom_Gbest[0] = algorythm1.nom_GB;
				for (i1 = 0; i1 < N[0]; i1++) {
					if (i1 < algorythm1.Nh)
						nom_x[0][i1] = algorythm1.nom_swarmH[i1];
					else
						nom_x[0][i1] = algorythm1.nom_swarmP[i1];
				}
		/*	}

#pragma omp section 
			{*/
			algorythm2.N = N[1];
			if (t2 >= 1)
			for (i2 = 0; i2 < N[1]; i2++) {
				for (d2 = 0; d2 < D; d2++) {
					algorythm2.x[i2][d2] = x[1][i2][d2];
				}
			}
			algorythm2.Ibav(ann);
			for (d2 = 0; d2 < D; d2++) {
				for (i2 = 0; i2 < N[1]; i2++)
					x[1][i2][d2] = algorythm2.x[i2][d2];
				Gbest[1][d2] = algorythm2.GBest[d2];
			}
			nom_Gbest[1] = algorythm2.nom_GBest;
			for (i2 = 0; i2 < N[1]; i2++)
				nom_x[1][i2] = algorythm2.nom_x[i2];
		/*}

			#pragma omp section
					{*/
					algorythm3.N = N[2];
					if (t >= 1)
					for (i = 0; i < N[2]; i++) {
						for (d = 0; d < D; d++) {
							algorythm2.x[i][d] = x[2][i][d];
						}
					}
					algorythm3.NSRaFAv(ann);
					for (d = 0; d < D; d++) {
						for (i = 0; i < N[2]; i++)
							x[2][i][d] = algorythm3.x[i][d];
						Gbest[2][d] = algorythm3.gBest[d];
					}
					nom_Gbest[2] = algorythm3.nom_gBest;
					for (i = 0; i < N[2]; i++)
						nom_x[2][i] = algorythm3.nom_x[i];
					/*}
					#pragma omp section
					{*/
					algorythm4.NP = N[3];
					if (t >= 1)
					for (i = 0; i < N[3]; i++) {
						for (d = 0; d < D; d++) {
							algorythm2.x[i][d] = x[3][i][d];
						}
					}
					algorythm4.CSKHv(ann);
					for (d = 0; d < D; d++) {
						for (i = 0; i < N[3]; i++)
							x[3][i][d] = algorythm4.x[i][d];
						Gbest[3][d] = algorythm4.gBest[d];
					}
					nom_Gbest[3] = algorythm4.nom_gBest;
					for (i = 0; i < N[3]; i++)
					nom_x[3][i] = algorythm4.nom_x[i];
					//}
		//}
	//}
		GLOBAL_COUNTER = algorythm1.GLOBAL_COUNTER + algorythm2.GLOBAL_COUNTER + algorythm3.GLOBAL_COUNTER + algorythm4.GLOBAL_COUNTER;
		
		/*for (j = 0; j < M; j++)
		{
			for (d = 0; d < D; d++)
				fout << Gbest[j][d] << '\t';
			fout << '\n';
		}*/
}

void COBRA::best_fitness()
{
	int i, j;
	for (j = 0; j < M; j++)
		average_fitness[j] = 0;

	//просчет средних пригодностей каждого алгоритма
	for (j = 0; j < M; j++)
	{
		for ( i = 0; i < N[j]; i++)
			average_fitness[j] += nom_x[j][i];
		average_fitness[j] = average_fitness[j] / (double) N[j];
	}

	best_algorythm = 0;
	for ( j = 1; j < M; j++)
		if (average_fitness[best_algorythm]>average_fitness[j]) best_algorythm = j;
}

void COBRA::population_size_change(ANN& ann)
{
	int i, j, d, sum = 0;
	for (j = 0; j < M; j++)
	if (j != best_algorythm && N[j] * 0.9 >= N_MIN)
	{
		sum += 0.1*N[j];
		N[j] = 0.9*N[j];
	}
	if (sum == 0) sum = M - 1;
	if ((N[best_algorythm] + sum) <= N_MAX)
		N[best_algorythm] = N[best_algorythm] + sum;
}

void COBRA::best()
{
	int i, d,k=0;
	for ( i = 0; i < M; i++)
	{
		if (nom_Gbest[i] < nom_GBest)
		{
			for (d = 0; d < D;d++)
				GBest[d] = Gbest[i][d];
			nom_GBest = nom_Gbest[i];
			k++;
		}
	}
	if (k == 0) time_without_сhange++;
}

void COBRA::exchange()
{
	int i, j, zxc;
	for ( j = 0; j < M; j++)
	{
		zxc = 0;
		for (i = N[j]-3; i <N[j]; i++)
		{
			if (j == zxc) zxc++;
			for (int d = 0; d < D;d++)
				x[j][i][d] = Gbest[zxc][d];
			nom_x[j][i] = nom_Gbest[zxc];
			zxc++;
		}
	}
	
}

void COBRA::COBRA_v(PSOPB_v & algorythm1, iba_v & algorythm2, NSRaFA_v& algorythm3, CSKH_v& algorythm4, ANN& ann)
{
	int j, i, d;
	int tt, TT;
	ofstream fout;
	//fout.open("QWE.txt");
	//int start_time = clock();

	GLOBAL_COUNTER = 0;

	time_without_сhange = 0;
	input(algorythm1, algorythm2, algorythm3, algorythm4, ann);
	TT = T;
	//d = omp_get_num_procs();
	//omp_set_num_threads(omp_get_num_procs()); //устанавливаем количество потоков равным количеству доступных процессоров в системе
//#pragma omp for ordered schedule(dynamic)
	for (tt = 0; tt < TT;tt++)
			{
				t = tt;
				algorythm1.t = t;
				algorythm2.t = t;
				algorythm3.t = t;
				algorythm4.t = t;

				//cout << "Итерация:	" << t << '\t' << "Время выполнения программы:	" << (double)(clock() - start_time) / 60000. << "мин" << '\t' << "Minimum:" << '\t' << nom_GBest << '\n';

				work_of_algorythms(algorythm1, algorythm2, algorythm3, algorythm4, ann);
				population_size_change(ann);
				best_fitness();
				best();

				if (GLOBAL_COUNTER > 100) {
					exchange();
				}
				if (t > 10)
				if (time_without_сhange != 0)
				for (j = 0; j < M; j++)
				{
					if (N[j] + K <= N_MAX)
					{
						N[j] += K;
						for (i = N[j] - K; i < N[j]; i++) {
							for (d = 0; d < D; d++)
								x[j][i][d] = ((double)rand() * (G[1][d] - G[0][d]) / RAND_MAX + G[0][d]); //[a; b]
							nom_x[j][i] = f(GLOBAL_COUNTER, x[j][i], ann);
						}
					}
				}
				else
				for (j = 0; j < M; j++)
				if (N[j] - K >= N_MIN)
					N[j] -= K;


				//cout << nom_GBest << '\n';


				/*fout2 << t << '\n';
				for (j = 0; j < M; j++)
				fout2 << N[j] << '\n';*/
			}
	
}