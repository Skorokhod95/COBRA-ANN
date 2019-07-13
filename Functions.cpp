#define _USE_MATH_DEFINES

#include<math.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>

#include "Functions.h"
#include "Algorithms.h"

//Функция пригодности
double f(int&C,double* x,ANN& ann)
{
	C++;
	ann.learning(x);
	return ann.Error;
}



//Функции для генерации случайности
double Cauchy()
{
	return tan(M_PI*((rand() / double(RAND_MAX)) - 0.5));
}
double Gaussian(double M, double D)
{
	/*Сложим n случайных чисел, используя стандартный ГСЧ :

	Согласно ЦПТ числа V образуют ряд значений, распределенный по нормальному закону.Эти числа тем лучше описывают нормальный закон,
	чем больше параметр n.На практике n берут равными 6 или 12. Заметим, что закон распределения чисел V имеет математическое
	ожидание mV = n / 2, σV = sqrt(n / 12).Поэтому он является смещенным относительно заданного произвольного.
	С помощью формулы z = (V – mV) / σV нормализуем этот ряд.Получим нормализованный закон нормального распределения чисел Z.То есть mz = 0, σz = 1.
	Формулой(сдвиг на mx и масштабирование на σx) преобразуем ряд Z в ряд x : x = z · σx + mx.*/

	double V = 0, z, mV, sigmaV;
	int i, n = 12;

	mV = n / 2.;
	sigmaV = sqrt(n / 12.);

	for (i = 0; i < n; i++)
		V += (rand() / double(RAND_MAX));

	z = (V - mV) / sigmaV;
	return z*D + M;
}
double Levy(double l)
{
	return pow((1. / (rand() / double(RAND_MAX))), (1. / l));
}


