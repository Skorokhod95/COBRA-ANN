#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include<omp.h>
#include <chrono>
#include "Algorithms.h"
#include "Functions.h"


#include "C:\Users\Shito\OneDrive\ДИПЛОМ\Бакалавр\Нейросеть\NeuralNet_class\Debug\Eigen\Dense"

using namespace Eigen;

using namespace std;
using namespace std::chrono;



void main()
{
	time_point<system_clock> start, end;
	start = steady_clock::now();
	setlocale(0, "");
	srand(time(0));
	ifstream fin("learn.txt"); // открыли файл для чтения
	ifstream test("test.txt"); // открыли файл для чтения
	ofstream fout3, fout44;
	fout3.open("out3.txt");
	
	COBRA cobra;
	PSOPB_v a1;
	iba_v a2;
	NSRaFA_v a3;
	CSKH_v a4;
	cobra.N_MAX = 400;
	cobra.N_MIN = 100;

	//int number_of_neurons, number_of_layers, number_of_outneurons, number_of_inputs;

	int i,j;
	
	ANN ann;
	ann.initialize();

	ann.number_of_layers = 1;
	ann.number_of_neurons = 5;
	ann.number_of_outneurons = 2;

	int m;
	char str[50];
	//В первой строке должно быть одно число, равное количеству примеров в файле (ну так как мне лень считать количество строк до конца файла)
	fin >> ann.max_learning_iteration;

	//Запоминаешь, где курсор
	//int gh = fin.tellg();

	//Считываешь первую строку полностью (размерность чартовской строки у меня наобум, так что советую побольше, если у тебя 60 показателей 
	//fin >> str;
	//Здесь он вроде как считает табуляцию что ли, в общем он у тебя выдает здесь количество элементов(уже целых чисел, а не как в чарте) в строке + 1
	//m = strlen(str)-1;
	ann.number_of_inputs = 57;

	//Здесь возвращаем курсор на место, чтобы считать заново "пробную" строку
	//fin.seekg(gh);
	//cout << ann.max_learning_iteration;
	ann.sample = new double*[ann.max_learning_iteration];
	for (i = 0; i < ann.max_learning_iteration; i++)
		ann.sample[i] = new double[ann.number_of_inputs+1];
	//Ну и считываем в массив все 
	for (i = 0; i < ann.max_learning_iteration; i++) {
		for (j = 0; j < ann.number_of_inputs + 1; j++) {
			fin >> ann.sample[i][j];
			//cout << ann.sample[i][j] << '\t';
		}
		//cout << '\n';
	}
		
	

	
	//РАССЧИТЫВАЕМ РАЗМЕРНОСТЬ ВЕКТОРА ДЛЯ ОПТИМИЗАЦИИ, ТЕ ПРОСЧИТЫВАЕМ КОЛИЧЕСТВО ВСЕХ СВЯЗЕЙ
	ann.D = (ann.number_of_layers - 1)*((int)pow((double)ann.number_of_neurons, 2.)) + ann.number_of_inputs*ann.number_of_neurons + ann.number_of_outneurons*ann.number_of_neurons + 1 - (1 + (ann.number_of_layers - 1)*ann.number_of_neurons);

	cobra.D = ann.D;

	
	ann.create_learning_sample();

	int MM;
	double **check, *n_class[2];

	test >> MM;

	check = new double*[MM];
	n_class[0] = new double[MM];
	n_class[1] = new double[MM];
	for (i = 0; i < MM; i++)
		check[i] = new double[ann.number_of_inputs];

	
	for (int g = 0; g < MM; g++) {
		for (int h = 0; h < ann.number_of_inputs; h++) {
			test >> check[g][h];
			check[g][h] = (check[g][h] - ann.x_min[h]) / (ann.x_max[h] - ann.x_min[h]);
		}
	}
	for (int g = 0; g < MM; g++)
		test >> n_class[0][g];

	ann.create_learning_sample();
	double Error;
	for (int k = 0; k < 5; k++) {

	cout << "Кобра №" << '\t' << MM<< '\n';
		
		cobra.COBRA_v(a1, a2, a3, a4, ann);

		for (int g = 0; g < MM; g++){
			n_class[1][g] = ann.work_of_ANN(cobra.GBest, check[g]);
		}
		Error = 0;
		for (int g = 0; g < MM; g++)
		{
			Error += fabs((n_class[1][g] - n_class[0][g]));
		}
		for (int g = 0; g < MM; g++) {
			//for (int h = 0; h < ann.number_of_inputs; h++)
				//fout3 << check[g][h] * (ann.x_max[h] - ann.x_min[h]) + ann.x_min[h] << '\t';
			fout3 << n_class[1][g] << '\n';
		}

		Error = Error / MM;
		cout << '\n' << (1 - Error) * 100 << '\n';
	}
	
	fout3.close();
	

	//////////////////////////////////////////////////
	
	end = steady_clock::now();
	cout << "Время выполнения: " << duration_cast<milliseconds>	(end - start).count() << "ms\n";
	system("pause");
}
