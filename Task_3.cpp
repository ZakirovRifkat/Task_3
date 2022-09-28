#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

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

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Задание 3. Проблема собственных значений. Вариант 4. Матрица - 8. \n";
	int n;
	double e = 0.000001, max = 0;
	cout << "Введите размерность матрицы А:\n"
		<<"n = ";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n));
	vector<double> exSolution, my_vector_1(n);

	double exact_number = 11;
	cout << "Заполните вашу матрицу:\n";
	for(int i=0; i<n; i++)
		for (int j = 0; j < n; j++)
			cin >> matrix_A[i][j];

	cout << "Исходная матрица:\n";
	show(matrix_A);

	//max = degree_method(matrix_A, exact_number);
	//my_vector_1 = own_vector(matrix_A, max);
	//cout << "Степенной метод макс.собст.число = " << max << "\n";
	////	<< "Соответствующий ему собственный вектор ( ";
	////for (int i = 0; i < n; i++)
	////{
	////	if (i != n-1)
	////		cout << my_vector_1[i] << " ; ";
	////	else
	////		cout << my_vector_1[i] << " )\n";
	////}
	return(0);
}