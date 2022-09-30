#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

using namespace std;

void show(vector<vector<double>> matrix)
{
	int n = matrix.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "|" << matrix[i][j];
		cout << "|\n";
	}
}
vector<double> multiplyMatrixVector(vector<vector<double>>A, vector<double>b)
{
	int n = A.size();
	vector<double>Ab(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Ab[i] += A[i][j] * b[j];
	return Ab;
}
vector<vector<double>> multiplyMatrixMatrix(vector <vector <double>> A, vector <vector <double>> B)
{
	int n = A.size();
	vector <vector <double>>AB(n, vector<double>(n));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				AB[i][j] += A[i][k] * B[k][j];
	return AB;
}
double scalar(vector<double> a, vector<double> b)
{
	double ab = 0;
	for (int i = 0; i < a.size(); i++)
		ab += a[i] * b[i];
	return ab;
}
double aposteriori_estimate(vector<vector<double>> A, vector<double> Y, double lambda)
{
	int n = Y.size();
	vector<double> ay(n);
	ay = multiplyMatrixVector(A, Y);
	for (int i = 0; i < n; i++)
		Y[i] = lambda * Y[i];
	for (int i = 0; i < n; i++)
		ay[i] = ay[i] - Y[i];
	return sqrt(scalar(ay,ay)) / sqrt(scalar(Y,Y));
}
vector<double> degree_method(vector<vector<double>> A, double e)
{
	int n = A.size(), step =0;
	double estimate = 1, lambda=0;
	vector<double> y_0(n), y_k(n), solution(n+2);
	y_0[0] = 1;
	for (int i = 1; i < n; i++)
		y_0[i] = 1;
	while (estimate > e)
	{
		y_k = multiplyMatrixVector(A, y_0);
		lambda = y_k[0] / y_0[0];
	/*	cout << "\nлямбда = " << lambda << "\n";
		cout << "\nВектор y_k:\n";
		for (int i = 0; i < n; i++)
			cout << y_k[i]<<"\n";*/
		y_0 = y_k;
		estimate = aposteriori_estimate(A, y_k, lambda);
		//cout << "estimate = " << estimate<<"\n";
		step++;
	}
	solution[0] = step;
	solution[1] = lambda;
	for (int i = 2; i < n+2; i++)
		solution[i] = y_k[i - 2];
	return solution;
}
vector<double> scalar_method(vector<vector<double>> A, double e)
{
	int n = A.size(), step = 0;
	double estimate = 1, lambda = 0;
	vector<double> y_0(n), y_k(n), solution(n + 2);
	y_0[0] = 1;
	for (int i = 1; i < n; i++)
		y_0[i] = 1;
	while (estimate > e)
	{
		y_k = multiplyMatrixVector(A, y_0);
		lambda = scalar(y_k, y_0) / scalar(y_0, y_0);
		/*	cout << "\nлямбда = " << lambda << "\n";
		cout << "\nВектор y_k:\n";
		for (int i = 0; i < n; i++)
			cout << y_k[i]<<"\n";*/
		y_0 = y_k;
		estimate = aposteriori_estimate(A, y_k, lambda);
		//cout << "estimate = " << estimate<<"\n";
		step++;
	}
	solution[0] = step;
	solution[1] = lambda;
	for (int i = 2; i < n + 2; i++)
		solution[i] = y_k[i - 2];
	return solution;
}
vector<double> reverse_scalar(vector<vector<double>> A, double lambda)
{
	int n = A.size();
	double reverse_lambda = 0;
	vector<vector<double>> B (n, vector<double>(n));
	vector<double> vector(n+2);
	for (int i = 0; i < n; i++)
		B[i][i] = A[i][i] - lambda;
	vector = scalar_method(B, 0.000001);
	reverse_lambda = vector[1] + lambda;
	vector[1] = reverse_lambda;
	return vector;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Задание 3. Проблема собственных значений. Вариант 4. Матрица 8.\n";
	int n;
	double e = 0.000001, ex_max_value = 1.110;
	ifstream file("\Text.txt");
	file >> n;
	cout << "Размерность матрицы n = " << n<<"\n";
	vector<vector<double>> matrix_A(n, vector<double>(n));
	vector<double> sol_1(n + 2), sol_2(n+2), sol_3(n+2),x(n);

	//Заполняем матрицу
	for(int i=0; i<n; i++)
		for (int j = 0; j < n; j++)
			file >> matrix_A[i][j];

	cout << "Исходная матрица:\n";
	show(matrix_A);

	//степенной метод
	sol_1 = degree_method(matrix_A, 0.001);
	cout << "\nМаксимальное собственное число cтепенным методом = " << (sol_1[1]) << "\n"
		<< "Соотвествующий ему собственный вектор: ( ";
	for (int i = 2; i < n + 2; i++)
	{
		if (i != n + 1)
			cout << sol_1[i] << "; ";
		else
			cout << sol_1[i] << " )";
	}
	cout << "; Кол-во итераций = " << sol_1[0]<<"\n";
;
	
	//метод скалярных произведений
	sol_2 = scalar_method(matrix_A, 0.000001);
	cout << "\nМаксимальное собственное число методом скалярных произведений = " << (sol_2[1]) << "\n"
		<< "Соотвествующий ему собственный вектор: ( ";
	for (int i = 2; i < n + 2; i++)
	{
		if (i != n + 1)
			cout << sol_2[i] << "; ";
		else
			cout << sol_2[i] << " )";
	}
	cout << "; Кол-во итераций = " << sol_2[0] << "\n";


	//обратная граница:
		sol_3 = reverse_scalar(matrix_A, sol_2[1]);
	cout << "\nПротивоположная граница спектра = " << (sol_3[1]) << "\n"
		<< "Соотвествующий ему собственный вектор: ( ";
	for (int i = 2; i < n + 2; i++)
	{
		if (i != n + 1)
			cout << sol_3[i] << "; ";
		else
			cout << sol_3[i] << " )";
	}
	cout << "; Кол-во итераций = " << sol_3[0] << "\n";

	return(0);
}