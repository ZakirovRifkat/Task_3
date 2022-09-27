#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

void show(vector<vector<double>> matrix)
{
	int n = matrix.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "\t" << matrix[i][j];
		cout << "\n";
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
double spectr_radius(vector<vector<double>> A)
{
	int n = A.size();
	double lambda = 0;
	vector<double> x_0(n), x_k(n);
	for (int i = 0; i < n; i++)
		x_0[i] = 1;
	for (int k = 0; k < 50; k++)
	{
		x_k = multiplyMatrixVector(A, x_0);
		lambda = scalar(x_k, x_0) / scalar(x_0, x_0);
		x_0 = x_k;
	}
	return fabs(lambda);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Задание 3. Проблема собственных значений. Вариант 4. Матрица - 8. \n";
	int n;
	cout << "Введите размерность матрицы А\n";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n));
	cout << "Заполните вашу матрицу:";
	for(int i=0; i<n; i++)
		for (int j = 0; j < n; j++)
			cin >> matrix_A[i][j];

	cout << "Исходная матрица:\n";
	show(matrix_A);
	return(0);
}