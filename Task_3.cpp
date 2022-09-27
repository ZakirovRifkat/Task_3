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
double lenght(vector<double> vector)
{
	double sum = 0;
	for (int i = 0; i < vector.size(); i++)
		sum += pow(vector[i], 2);
	return sqrt(sum);
}

//vector<double> Gauss(vector<vector<double>> matrix_A, vector<double> vector_b)
//{
//	int p, n = matrix_A.size();
//	double r, c, s;
//	vector<double>x(n), b(n);
//	vector<vector<double>>a(n, vector<double>(n));
//	a = matrix_A;
//	b = vector_b;
//	for (int k = 0; k < n; k++)
//	{
//		p = k;
//		for (int m = k + 1; m < n; m++)
//			if (abs(a[p][k]) < abs(a[m][k])) //поиск максимального ведущего элемента
//				p = m;
//		for (int j = k; j < n; j++)
//		{
//			r = a[k][j];
//			a[k][j] = a[p][j];   //перестановка строк
//			a[p][j] = r;
//		}
//		r = b[k];
//		b[k] = b[p];   //перестановка свободных членов
//		b[p] = r;
//		for (int m = k + 1; m < n; m++)
//		{
//			c = a[m][k] / a[k][k];
//			b[m] = b[m] - c * b[k]; //приведение матрицы к верхнетреугольному виду
//			for (int i = k; i < n; i++)
//				a[m][i] = a[m][i] - c * a[k][i];
//		}
//	}
//	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
//	for (int k = n - 1; k >= 0; k--)
//	{
//		s = 0;
//		for (int i = k + 1; i < n; i++)				//обратный ход метода Гаусса
//			s = s + a[k][i] * x[i];
//		x[k] = (b[k] - s) / a[k][k];
//	}
//	return x;
//}

double scalar(vector<double> a, vector<double> b)
{
	double ab = 0;
	for (int i = 0; i < a.size(); i++)
		ab += a[i] * b[i];
	return ab;
}
double degree_method(vector<vector<double>> A, double exSol)
{
	int n = A.size();
	double lambda = 0;
	vector<double> x_0(n), x_k(n);
	for (int i = 0; i < n; i++)
		x_0[i] = 1;
	/*while(fabs(lambda - exSol) > 0.001)*/
	for(int i = 0; i<10; i++)
	{
		x_k = multiplyMatrixVector(A, x_0);
		lambda = scalar(x_k, x_0) / scalar(x_0, x_0);
		x_0 = x_k;
	}
	return lambda;
}
vector<double> own_vector(vector<vector<double>> matrix, double lambda)
{
	int n = matrix.size();
	double l = 0;
	vector<vector<double>> A(n, vector<double>(n));
	vector<double> null(n), own(n);
	for (int i = 0; i < n; i++)
		null[i] = 0;
	A = matrix;
	for (int i = 0; i < n; i++)
		A[i][i] = A[i][i] - lambda;
	//own = Gauss(A,null); ???????????????????????????????????????
	l = lenght(own);
	if (l != 1)
		for (int i = 0; i < n; i++)
			own[i] = own[i] / l;
	return own;
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

	max = degree_method(matrix_A, exact_number);
	my_vector_1 = own_vector(matrix_A, max);
	cout << "Степенной метод макс.собст.число = " << max << "\n";
	//	<< "Соответствующий ему собственный вектор ( ";
	//for (int i = 0; i < n; i++)
	//{
	//	if (i != n-1)
	//		cout << my_vector_1[i] << " ; ";
	//	else
	//		cout << my_vector_1[i] << " )\n";
	//}
	return(0);
}