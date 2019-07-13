#define _USE_MATH_DEFINES

#include "Functions.h"
#include "Algorithms.h"
#include <fstream>
#include <iostream>
#include "C:\Users\Shito\OneDrive\ÄÈÏËÎÌ\Áàêàëàâğ\Íåéğîñåòü\NeuralNet_class\Debug\Eigen\Dense"

using namespace std;
using namespace Eigen;

void ANN::initialize()
{
	weight_out = new double[D];

	sample = new double*[max_learning_iteration];
	for (i = 0; i < max_learning_iteration; i++)
		sample[i] = new double[number_of_inputs + number_of_outneurons];
	int i;
	for (i = 0; i < D; i++)
	{
		weight_out[i] = ((double)rand() * (b - a) / RAND_MAX + a);
	}
}

void ANN::create_learning_sample()
{
	int i, j, k;

	x_max = new double[number_of_outneurons+number_of_inputs];
	x_min = new double[number_of_outneurons+number_of_inputs];

	/*MatrixXd training(max_learning_iteration, number_of_inputs + number_of_outneurons);

	for (i = 0; i < max_learning_iteration; i++)
	{
		for (j = 0; j < number_of_inputs; j++)
			training(i, j) = ((double)rand() * (b - a) / RAND_MAX + a);
		for (k = number_of_inputs; k < number_of_inputs + number_of_outneurons; k++) {
			training(i, k) = training_function(training.row(i).segment(0, number_of_inputs), number_of_inputs);
			if (i == 0) {
				x_max[k - number_of_inputs] = training(i, k);
				x_min[k - number_of_inputs] = training(i, k);
			}
			else {
				if (x_max[k - number_of_inputs] < training(i, k)) x_max[k - number_of_inputs] = training(i, k);
				if (x_min[k - number_of_inputs] > training(i, k)) x_min[k - number_of_inputs] = training(i, k);
			}
		}
	}

	for (i = 0; i < max_learning_iteration; i++) {
		for (j = 0; j < number_of_inputs; j++)
			sample[i][j] = training(i, j);
		for (j = number_of_inputs; j < number_of_inputs + number_of_outneurons; j++)
			sample[i][j] = (training(i, j)-x_min[j-number_of_inputs]) / (x_max[j-number_of_inputs]-x_min[j-number_of_inputs]);
	}*/

	for (i = 0; i < max_learning_iteration; i++)
	{
		for (k = 0; k < number_of_inputs + 1; k++) {
			if (i == 0) {
				x_max[k] = sample[i][k];
				x_min[k] = sample[i][k];
			}
			else {
				if (x_max[k] < sample[i][k]) 
					x_max[k] = sample[i][k];
				if (x_min[k] > sample[i][k]) 
					x_min[k] = sample[i][k];
			}
		}
	}
	//for (k = 0; k < number_of_inputs + number_of_outneurons; k++)
		//cout << x_min[k] << '\t';
	//cout << '\n';

	//ÍÎĞÌÈĞÎÂÊÀ
	for (i = 0; i < max_learning_iteration; i++) {
		for (j = 0; j < number_of_inputs + number_of_outneurons; j++)
		{
			sample[i][j] = (sample[i][j] - x_min[j]) / (x_max[j] - x_min[j]);
		}
	}
}

double ANN::activation_function(int variant_of_function, double x)
{
	int beta = 1;

	switch (variant_of_function)
	{
	case 1:	//ñèãìîèä
	{			
		if (x<-10) return 0;
		else if (x>10) return 1;
		else return 1. / (1 + exp(-beta*x));
	}
	case 2:	//áèïîëÿğíàÿ ôóíêöèÿ
		return	tanh(beta*x);
	case 3:
	{
			  if (x >= 0) return 1;
			  else return 0;
	}
	}
}


double ANN::training_function(VectorXd x, int D)
{
	int i, j, d;
	double sum = 0;
	for (j = 0; j < D; j++)
		sum += pow(x(j), 2.);
	return sum;

	/*sum = 0;
	for (d = 0; d < D; d++)
		sum += (pow(x[d], 2.) - 10 * cos(2 * M_PI*x[d]));
	return 200 + sum;*/

}

void ANN::learning(double *weight_out)
{
	double Error1=0;
	Error = 0;

	//ofstream output;
	//output.open("error.txt");

	/*MatrixXd new_matrix(êîëè÷åñòâî ñòğîê, êîëè÷åñòâî ñòîëáöîâ);*/
	MatrixXd weight_all(number_of_neurons, number_of_neurons*(number_of_layers - 1)), weight(number_of_neurons, number_of_neurons);
	MatrixXd input_weight(number_of_inputs, number_of_neurons), output_weight(number_of_neurons, number_of_outneurons);
	VectorXd neurous_out(number_of_neurons), neurous_in(number_of_neurons), out(number_of_outneurons), in(number_of_inputs);
	MatrixXd training(max_learning_iteration, number_of_inputs + number_of_outneurons);



	//ÇÀÏÈÑÜ ÂÅÊÒÎĞÀ ÂÅÑÎÂ Â ÌÀÒĞÈÖÓ
	int black = 0, num_w = 0;
	//ÇÀÏÈÑÜ ÂÅÊÒÎĞÀ ÂÅÑÎÂ Â ÌÀÒĞÈÖÓ
	//Çàïèñü âåñîâ íà÷àëüíîãî ñëîÿ
	for (i = 0; i < number_of_inputs; i++) {
		for (j = 0; j < number_of_neurons; j++) {
			input_weight(i, j) = weight_out[num_w];
			if (j % ((number_of_neurons - 1) + black*number_of_neurons) == 0 && j != 0) {
				input_weight(i, j) = 0;
				black++;
				num_w--;
			}
			num_w++;
		}
		black = 0;
	}
	black = 0;
	//Çàïèñü âåñîâ âíóòğåííèõ ñëîåâ
	for (i = 0; i < number_of_neurons; i++) {
		for (j = 0; j < number_of_neurons*(number_of_layers - 1); j++) {
			weight_all(i, j) = weight_out[num_w];
			if (j % ((number_of_neurons - 1) + black*number_of_neurons) == 0 && j != 0) {
				weight_all(i, j) = 0;
				black++;
				num_w--;
			}
			num_w++;
		}
		black = 0;
	}
	//Çàïèñü âåñîâ ïîñëåäíåãî ñëîÿ
	for (i = 0; i < number_of_neurons; i++) {
		for (j = 0; j < number_of_outneurons; j++) {
			output_weight(i, j) = weight_out[num_w];
			num_w++;
		}
	}


	double aaaa;
	for (i = 0; i < max_learning_iteration; i++) {
		for (j = 0; j < number_of_inputs + 1; j++) {
			training(i, j) = sample[i][j];
			}
	}

	for (learning_iteration = 0; learning_iteration < max_learning_iteration; learning_iteration++) {
		
		//ÍÀ×ÀËÎ ĞÀÁÎÒÛ (ÏÅĞÂÛÉ ÑËÎÉ)
		in.transpose() = training.row(learning_iteration).segment(0, number_of_inputs).transpose();
		
		neurous_in.transpose() = in.transpose()*input_weight;
		for (i = 0; i < number_of_neurons - 1; i++) {
			neurous_out(i) = activation_function(1, neurous_in(i));
		}
		neurous_out(number_of_neurons - 1) = 1;

		//ĞÀÁÎÒÀ ÂÍÓÒĞÈ ÑÅÒÈ
		for (j = 0; j < number_of_layers - 1; j++) {
			for (i = 0; i < number_of_neurons; i++)
				weight.col(i) = weight_all.col(i + number_of_neurons*j);
			//neurous_out()
			neurous_in.transpose() = neurous_out.transpose()*weight;
			for (i = 0; i < number_of_neurons - 1; i++) {
				neurous_out(i) = activation_function(1, neurous_in(i));
			}
			neurous_out(number_of_neurons - 1) = 1;
		}

		//ĞÀÁÎÒÀ ÏÎÑËÅÄÍÅÃÎ ÍÅÉĞÎÍÀ
		out.transpose() = neurous_out.transpose()*output_weight;
		for (i = 0; i < number_of_outneurons; i++)
			out(i) = activation_function(1, out(i));

		//output << training(learning_iteration, 0) << '\t' << out << '\t' << training(learning_iteration, 1) << '\n';
		/*for (i = 0; i < number_of_outneurons; i++)
			Error += pow(training(learning_iteration, number_of_inputs + i) - out(i), 2.);*/
		double maxim1, maxim2=0;
		maxim1 = out(0);
		for (i = 0; i < number_of_outneurons; i++)
		if (maxim1>out(i)){
			maxim1 = out(i);
			maxim2 = i;
		}

		/*for (i = 0; i < number_of_outneurons; i++) {
			if (out(i)>0.5) out(i) = 1;
			else out(i) = 0;*/
			Error += fabs(training(learning_iteration, number_of_inputs) - maxim2);
		//}
	}

	

	Error = Error/max_learning_iteration * 100;

}

double ANN::work_of_ANN(double *weight_out, double* input)
{
	int i, j;

		
	/*MatrixXd new_matrix(êîëè÷åñòâî ñòğîê, êîëè÷åñòâî ñòîëáöîâ);*/
	MatrixXd weight_all(number_of_neurons, number_of_neurons*(number_of_layers - 1)), weight(number_of_neurons, number_of_neurons);
	MatrixXd input_weight(number_of_inputs, number_of_neurons), output_weight(number_of_neurons, number_of_outneurons);
	VectorXd neurous_out(number_of_neurons), neurous_in(number_of_neurons), out(number_of_outneurons), neurous_out_one(number_of_inputs);
	MatrixXd training(max_learning_iteration, number_of_inputs + number_of_outneurons);

	//ÇÀÏÈÑÜ ÂÅÊÒÎĞÀ ÂÅÑÎÂ Â ÌÀÒĞÈÖÓ
	int black = 0, num_w = 0;
	//ÇÀÏÈÑÜ ÂÅÊÒÎĞÀ ÂÅÑÎÂ Â ÌÀÒĞÈÖÓ
	//Çàïèñü âåñîâ íà÷àëüíîãî ñëîÿ
	for (i = 0; i < number_of_inputs; i++) {
		for (j = 0; j < number_of_neurons; j++) {
			input_weight(i, j) = weight_out[num_w];
			if (j % ((number_of_neurons - 1) + black*number_of_neurons) == 0 && j != 0) {
				input_weight(i, j) = 0;
				black++;
				num_w--;
			}
			num_w++;
		}
		black = 0;
	}
	black = 0;
	//Çàïèñü âåñîâ âíóòğåííèõ ñëîåâ
	for (i = 0; i < number_of_neurons; i++) {
		for (j = 0; j < number_of_neurons*(number_of_layers - 1); j++) {
			weight_all(i, j) = weight_out[num_w];
			if (j % ((number_of_neurons - 1) + black*number_of_neurons) == 0 && j != 0) {
				weight_all(i, j) = 0;
				black++;
				num_w--;
			}
			num_w++;
		}
		black = 0;
	}
	//Çàïèñü âåñîâ ïîñëåäíåãî ñëîÿ
	for (i = 0; i < number_of_neurons; i++) {
		for (j = 0; j < number_of_outneurons; j++) {
			output_weight(i, j) = weight_out[num_w];
			num_w++;
		}
	}

	//ÍÀ×ÀËÎ ĞÀÁÎÒÛ (ÏÅĞÂÛÉ ÑËÎÉ)
	
	for (i = 0; i < number_of_inputs; i++)
		neurous_out_one(i) = input[i];
	neurous_in.transpose() = neurous_out_one.transpose()*input_weight;
	for (i = 0; i < number_of_neurons - 1; i++) {
		neurous_out(i) = activation_function(1, neurous_in(i));
	}
	neurous_out(number_of_neurons - 1) = 1;

	//ĞÀÁÎÒÀ ÂÍÓÒĞÈ ÑÅÒÈ
	for (j = 0; j < number_of_layers - 1; j++) {
		for (i = 0; i < number_of_neurons; i++)
			weight.col(i) = weight_all.col(i + number_of_neurons*j);
		neurous_in.transpose() = neurous_out.transpose()*weight;
		for (i = 0; i < number_of_neurons - 1; i++) {
			neurous_out(i) = activation_function(1, neurous_in(i));
		}
		neurous_out(number_of_neurons - 1) = 1;
	}

	//ĞÀÁÎÒÀ ÏÎÑËÅÄÍÅÃÎ ÍÅÉĞÎÍÀ
	out.transpose() = neurous_out.transpose()*output_weight;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	for (i = 0; i < number_of_outneurons; i++)
		out(i) = activation_function(1, out(i));

	double maxim1, maxim2=0;
	maxim1 = out(0);
	for (i = 0; i < number_of_outneurons; i++)
	if (maxim1>out(i)){
		maxim1 = out(i);
		maxim2 = i;
	}

	return maxim2;
}