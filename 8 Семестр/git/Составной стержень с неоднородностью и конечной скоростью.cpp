#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <omp.h>

#include <shlwapi.h>
#include <windows.h>
#pragma comment (lib, "shell32.lib ")

using namespace std;

int type = 0;

void printm(double*b, int n); // Вывод вектора
void printm(double**a, int n);
double *m3mult(double *a1, double *a2, double *a3, double *b, int n); // Произведение трехдиагональной "матрицы" на вектор
double *mmult(double a, double *b, int n); // Произведение константы на вектор
double *mmult(double **a, double *b, double *c, int n); // Для критерия остановы (невязка ||Ax-b|| < eps)
double *sum(double *a, double *b, int n); // Сумма векторов
void equate(double *a, double *b, int n); // Приравнивание векторов
double q1(double t); // Поток слева
double q2(double t); // Поток справа
double func(double x, double t); // Функция для теста
double ts(double t, double alpha); // Температура среды
double kyb(double *x, int n); // Норма вектора кубическая

void make_P(double *P, double tau, double t_relax, double **K, double *T_new, double *T_old, int n); // Реккурентная формула для вектора P
void make_rhs(double *d, double tau, double t_relax, double *f, double *P, double **rhs_dop, double *T, int n); // Правая часть для счета
void mmult_nloc(double *b1, double **a, double *b, int n); // Произведение матрицы вида матрицы нелокальности (с известным числом диагоналей) на вектор
void zeidel(double**a, double **aij_ii, double*b, double*x, int n); // Метод Зейделя. aij_ii[i][j] = a[i][j] / a[i][i]


//Радиус нелокальности (в элементах)
int k = 3;

int main()
{
	setlocale(LC_ALL, "Russian");
	
	double L = 1.; // Параметры стержня

	//double ro_l = 5400.; double c_l = 600.; // Нитрид титана
	//double la_l = 41.8; // Коэффициент теплопроводности на левой части стержня

	//double ro_r = 8940.; double c_r = 485.; // Медь
	//double la_r = 401; // Коэффициент теплопроводности на правой части стержня

	//double ro_r = 7800.; double c_r = 460.; // Сталь
	//double la_r = 22.4; // коэффициент теплопроводности на правой части стержня
		
	double ro_l = 1.; double c_l = 1.;
	double la_l = 1.; // Коэффициент теплопроводности на левой части стержня
	double ro_r = 3.; double c_r = 3.;
	double la_r = 3.; // Коэффициент теплопроводности на правой части стержня

	int n = 100; // Кол-во конечных элементов
	double h = L / n; // Шаг по стержню
	n++; // Количество узлов на единицу больше количества элементов	
	double tau = 0.01; // Шаг по времени
	double t1 = 4; // Время
	double m = round(t1 / tau); // Количество шагов по времени
	double T0 = 400; // Начальное условие
	int n1 = 5 * n / 10; // Точка соединения двух частей стержня
	int n2 = n - n1 + 1;
	//cout << n1 << " " << n2 << endl;
	
	double t_relax = 1;

	double p1 = 0.5;
	double p2 = 1 - p1;  // Параметры нелокальности для левой части стержня
	double a = (double) k;

	double *c1 = new double[n];
	double *c2 = new double[n];		
	double *c3 = new double[n];
								// Создаем вектора трёх диагоналей матриц С и К 
	double *k1 = new double[n];
	double *k2 = new double[n];	
	double *k3 = new double[n];

	// Считаем коэффициенты матриц для левой части стержня

	c1[0] = 0; c2[0] = ro_l*c_l*h / 3; c3[0] = ro_l*c_l*h / 6;
	k1[0] = 0; k2[0] = p1*la_l / h; k3[0] = -p1*la_l / h;

	for (int i = 1; i < n1 - 1; i++)
	{
		c1[i] = ro_l*c_l*h / 6;
		c2[i] = 2 * ro_l*c_l*h / 3;
		c3[i] = ro_l*c_l*h / 6; 				
		k1[i] = -p1*la_l / h;
		k2[i] = 2 * p1*la_l / h;
		k3[i] = -p1*la_l / h;
	}

	// В месте соединения
	c1[n1 - 1] = ro_l*c_l*h / 6; c2[n1 - 1] = ro_l*c_l*h / 3 + ro_r*c_r*h / 3; c3[n1 - 1] = ro_r*c_r*h / 6;;
	k1[n1 - 1] = -p1*la_l / h; k2[n1 - 1] = p1*la_l / h + la_r / h; k3[n1 - 1] = -la_r / h;

	// Считаем коэффициенты для правой части стержня
	for (int i = n1; i < n - 1; i++)
	{
		c1[i] = ro_r*c_r*h / 6;
		c2[i] = 2 * ro_r*c_r*h / 3;
		c3[i] = ro_r*c_r*h / 6; 				
		k1[i] = -la_r / h;
		k2[i] = 2*la_r / h;
		k3[i] = -la_r / h;
	}

	c1[n - 1] = ro_r*c_r*h / 6; c2[n - 1] = ro_r*c_r*h / 3; c3[n - 1] = 0;
	k1[n - 1] = -la_r / h; k2[n - 1] = la_r / h; k3[n - 1] = 0;

	// Запустить Mathematica, чтобы посчитать матрицу нелокальности
	/*cout << "n1 = " << n1 << ", a = " << a << endl;
	ShellExecuteA(NULL, "open", "\"A:\\Учеба\\Курсач\\8 Семестр\\C++\\Матрица.nb\"", NULL, NULL, SW_SHOW);
	system("pause");*/

	// Матрица нелокальности, которая впоследствии станет матрицей A для Зейделя

	double **matrixA = new double*[n];
	for (int i = 0; i < n; i++)
		matrixA[i] = new double[n];
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		matrixA[i][j] = 0;

	/*cout << "C:" << endl;
	printm(c3, n1);
	cout << endl;
	printm(c2, n1);
	cout << endl;
	printm(c1, n1);
	cout << endl;
	system("pause");*/

	ifstream input("A:/Учеба/Курсач/8 Семестр/C++/matrix_nloc.dat");
	if (input.is_open())
	for (int i = 0; i < n1; i++)
	for (int j = 0; j < n1; j++)
	{
		input >> matrixA[i][j];
		matrixA[i][j] *= la_l*p2 / (2 * a*h);
	}
	else
		cout << "Ошибка чтения файла" << endl;

	//printm(matrixA, n1);
	input.close();

	// Соединяем матрицы K
	matrixA[0][0] += k2[0];
	matrixA[0][1] += k3[0];
	for (int i = 1; i < n - 1; i++)
	{
		matrixA[i][i - 1] += k1[i];
		matrixA[i][i] += k2[i];
		matrixA[i][i + 1] += k3[i];
	}
	matrixA[n - 1][n - 1] += k2[n - 1];
	matrixA[n - 1][n - 2] += k1[n - 1];
		
	double **matrixK = new double*[n]; // Сохраним матрицу K для вычисления вектора P
	for (int i = 0; i < n; i++)
		matrixK[i] = new double[n];

	double dop = -t_relax*(1 - exp(-tau / t_relax)); // Коэффициент для матрицы K+K_nloc при сборке правой части для решения системы

	double **rhs_dop = new double*[n]; // Собираем C-(K+K_nloc)*dop (минус учтен в dop)
	for (int i = 0; i < n; i++)
		rhs_dop[i] = new double[n];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	{
		rhs_dop[i][j] = matrixA[i][j] * dop;
		matrixK[i][j] = matrixA[i][j];		// Сохраняем K
	}

	rhs_dop[0][0] += c2[0];
	rhs_dop[0][1] += c3[0];
	for (int i = 1; i < n - 1; i++)
	{
		rhs_dop[i][i - 1] += c1[i];
		rhs_dop[i][i] += c2[i];
		rhs_dop[i][i + 1] += c3[i];
	}
	rhs_dop[n - 1][n - 1] += c2[n - 1];
	rhs_dop[n - 1][n - 2] += c1[n - 1];


	double K_coef = tau - t_relax*(1 - exp(-tau / t_relax));



	// Собираем матрицу А для счета
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		matrixA[i][j] *= K_coef; // K_coef*K

	// K_coef*K + C
	matrixA[0][0] += c2[0];
	matrixA[0][1] += c3[0];
	for (int i = 1; i < n - 1; i++)
	{
		matrixA[i][i - 1] += c1[i];
		matrixA[i][i] += c2[i];
		matrixA[i][i + 1] += c3[i];
	}
	matrixA[n - 1][n - 1] += c2[n - 1];
	matrixA[n - 1][n - 2] += c1[n - 1];

	double **aij_ii = new double*[n];
	for (int i = 0; i < n; i++)
		aij_ii[i] = new double[n];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		aij_ii[i][j] = matrixA[i][j] / matrixA[i][i];

	double *T = new double[n]; // Вектор решения
	for (int i = 0; i < n; i++)
		T[i] = T0;

	cout << "\t\t Выберите тип задачи: " << endl;
	cout << "1. Постоянные потоки " << endl << "2. Граничные условия 2-го рода " << endl << "3. Граничные условия 2-го и 3-го рода" << endl << "Для выхода нажмите 0" << endl;
	cin >> type;
	if (type != 0) // Тело программы 
	{
		if (type == 1) // Постоянные потоки
		{
			double q1 = 100.; double q2 = 0; // Начальные условия
			double *f = new double[n]; // Вектор правой части
			f[0] = q1; f[n - 1] = q2;
			for (int i = 1; i < n - 1; i++)			// Заполняем вектор правой части
				f[i] = 0;			

			double *P = new double[n]; // Вспомогательный вектор, который считается по реккурентной формуле
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			double t1, t2;

			t1 = omp_get_wtime();

			ofstream fout;											// Для записи во временным слоям
			fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT14.dat");		

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				for (int i = 0; i < n; i++)
					fout << T[i] << " ";			// Для записи во временным слоям
				fout << endl;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);
			}
			t2 = omp_get_wtime();
			cout << "Время: " << t2 - t1 << " c" << endl;

			/*ofstream fout;
			fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT10.dat");*/
			//fout.open("G:/Курсач/C++/vecT1.dat"); 
			//fout.open("F:/Курсач/C++/vecT1.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] f;
			delete[] rhs;
			delete[] T_old;
			delete[] P;
		}

		if (type == 2) // Граничные условия 2-го рода
		{
			double *f = new double[n]; // Вектор правой части
			double t;

			double *P = new double[n]; // Вспомогательный вектор, который считается по реккурентной формуле
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			double t1, t2;

			t1 = omp_get_wtime();

			//ofstream fout;											// Для записи во временным слоям
			//fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT2.dat");

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 - q2(t);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);

				//for (int i = 0; i < n; i++)
				//	fout << T[i] << " ";			// Для записи во временным слоям
				//fout << endl;
			}

			t2 = omp_get_wtime();
			cout << "Время: " << t2 - t1 << " c" << endl;

			ofstream fout;
			fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT2.dat");
			//fout.open("F:/Курсач/C++/vecT2.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();
			
			delete[] f;
			delete[] rhs;
			delete[] T_old;
			delete[] P;
		}
		

		if (type == 3) // Граничные условия 2-го и 3-го рода
		{
			double *f = new double[n]; // Вектор правой части
			double t;

			double alpha = 10; // Коэффициент конвективного теплообмена

			matrixA[n - 1][n - 1] += tau*alpha;
			for (int i = 0; i < n; i++)
				aij_ii[i][n - 1] = matrixA[i][n - 1] / matrixA[i][i];
			for (int j = 0; j < n; j++)
				aij_ii[n - 1][j] = matrixA[n - 1][j] / matrixA[n - 1][n - 1];

			double *P = new double[n]; // Вспомогательный вектор, который считается по реккурентной формуле
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			//ofstream fout;											// Для записи во временным слоям
			//fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT3.dat");

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 + alpha*ts(t, alpha);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);

				//for (int i = 0; i < n; i++)
				//	fout << T[i] << " ";			// Для записи во временным слоям
				//fout << endl;
			}

			ofstream fout;
			fout.open("A:/Учеба/Курсач/8 Семестр/C++/vecT3.dat");
			//fout.open("F:/Курсач/C++/vecT3.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] f;
			delete[] rhs;
			delete[] T_old;
			delete[] P;
		}
		//printm(T, n);
	}

	system("pause");
	for (int i = 0; i < n; i++)
	{
		delete[] matrixA[i];
		delete[] aij_ii[i];
		delete[] rhs_dop[i];
		delete[] matrixK[i];
	}
	delete[] matrixA;
	delete[] aij_ii;
	delete[] rhs_dop;
	delete[] matrixK;

	delete[] c1;
	delete[] c2;
	delete[] c3;
	delete[] k1;
	delete[] k2;
	delete[] k3;
}//main

void zeidel(double**a, double **aij_ii, double*b, double*x, int n)
{
	double  norma_xr;

	double *y, *x0, *xr;
	double epsilon;

	epsilon = 1e-7;

	y = new double[n];
	x0 = new double[n];
	xr = new double[n];

	for (int i = 0; i < n; i++)
	{
		y[i] = b[i] / a[i][i];
		x0[i] = x[i];
	}

	int iter = 0;
	do
	{
		double dop;
		for (int i = 0; i < k + 1; i++)
		{
			x[i] = y[i];

			dop = 0;
			for (int j = i + 1; j < i + k + 2; j++)
				dop -= aij_ii[i][j] * x[j];

			for (int j = 0; j < i; j++)
				dop -= aij_ii[i][j] * x[j];

			x[i] += dop;

			xr[i] = fabs(x0[i] - x[i]);
			x0[i] = x[i];
		}


		for (int i = k + 1; i < n - k - 1; i++)
		{
			x[i] = y[i];

			dop = 0;
			for (int j = i + 1; j < i + k + 2; j++)
				dop -= aij_ii[i][j] * x[j];

			for (int j = i - k - 1; j < i; j++)
				dop -= aij_ii[i][j] * x[j];

			x[i] += dop;

			xr[i] = fabs(x0[i] - x[i]);
			x0[i] = x[i];
		}


		for (int i = n - k - 1; i < n; i++)
		{
			x[i] = y[i];

			dop = 0;
			for (int j = i + 1; j < n; j++)
				dop -= aij_ii[i][j] * x[j];

			for (int j = i - k - 1; j < i; j++)
				dop -= aij_ii[i][j] * x[j];

			x[i] += dop;

			xr[i] = fabs(x0[i] - x[i]);
			x0[i] = x[i];
		}

		iter++;

		norma_xr = kyb(xr, n);
		//} while (iter < 100);
	} while (norma_xr > epsilon);

	//cout << iter << endl;

	/*cout << "Решено за " << iter << " итераций " << endl;
	system("pause");*/

	delete[] y;
	delete[] x0;
	delete[] xr;
}

void mmult_nloc(double *b1, double **a, double *b, int n)
{
	for (int i = 0; i < k + 1; i++)
	{
		b1[i] = 0;

		for (int j = 0; j < i + k + 2; j++)

			b1[i] += a[i][j] * b[j];
	}

	for (int i = k + 1; i < n - k - 1; i++)
	{
		b1[i] = 0;

		for (int j = i - k - 1; j < i + k + 2; j++)

			b1[i] += a[i][j] * b[j];
	}

	for (int i = n - k - 1; i < n; i++)
	{
		b1[i] = 0;

		for (int j = i - k - 1; j < n; j++)

			b1[i] += a[i][j] * b[j];
	}
}

void make_rhs(double *d, double tau, double t_relax, double *f, double *P, double **rhs_dop, double *T, int n)
{
	double *rhs_dop_vec = new double[n];

	mmult_nloc(rhs_dop_vec, rhs_dop, T, n);

	double tau_t_relax = tau / t_relax;

	for (int i = 0; i < n; i++)
		d[i] = tau*(f[i] + exp(-tau_t_relax)*P[i]) + rhs_dop_vec[i];

	delete[] rhs_dop_vec;
}

void make_P(double *P, double tau, double t_relax, double **K, double *T_new, double *T_old, int n)
{
	double exp_dop = exp(-tau / t_relax);
	double dop = t_relax / tau*(1 - exp(-tau / t_relax));

	double *T = new double[n];
	for (int i = 0; i < n; i++)
		T[i] = dop*(T_new[i] - T_old[i]);

	double *T_dop = new double[n];
	mmult_nloc(T_dop, K, T, n);

	for (int i = 0; i < n; i++)
		P[i] = exp_dop*P[i] + T_dop[i];

	delete[] T, T_dop;
}

double *m3mult(double *a1, double *a2, double *a3, double *b, int n)
{
	double *b1;
	b1 = new double[n];
	b1[0] = a2[0] * b[0] + a3[0] * b[1];
	for (int i = 1; i < n - 1; i++)
		b1[i] = a1[i] * b[i - 1] + a2[i] * b[i] + a3[i] * b[i + 1];
	b1[n - 1] = a1[n - 1] * b[n - 2] + a2[n - 1] * b[n - 1];
	return b1;
	delete[] b1;
}

double *mmult(double a, double *b, int n)
{
	double *b1;
	b1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = a*b[i];
	}
	return b1;
	delete[] b1;
}

double *mmult_min(double **a, double *b, double *c, int n)
{
	double *b1;
	b1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = 0;

		for (int j = 0; j < n; j++)

			b1[i] += a[i][j] * b[j];

		b1[i] -= c[i];
	}
	return b1;
	delete[] b1;
}

double *sum(double *a, double *b, int n)
{
	double *b1;
	b1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = a[i] + b[i];
	}
	return b1;
	delete[] b1;
}

void printm(double*b, int n)
{
	for (int i = 0; i < n; i++)
		cout << b[i] << " ";
	cout << endl;
}

void printm(double**a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << a[i][j] << " ";
		cout << endl;
	}
}

void equate(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = b[i];
}

double q1(double t)
{
	if (type == 2)
		//return 1000;
		//return 4 * t*t*exp(-2.0*t);
		return t;
		//return pow(t, 5);
	if (type == 3)
		return t;
		/*return pow(t, 10);*/
		//return 100;
}

double q2(double t)
{
	//return 0;
	return t/exp(1.0);
}

//Правая часть
double func(double x, double t)
{
	if (type == 2)
		return exp(-x*x)*(1 - x + 2*t*(1 + x)*(1 + 2*x*(-2 + x)));
	//return 0;
	if (type == 3)
		return exp(-x*x)*(2 - x + t*(4 + 2 * x*(-3 + 2 * (-2 + x)*x)));
		//return 0; 
}

double ts(double t, double alpha)
{
	//return 0;
	return (t / exp(1.0))*(alpha - 3) / alpha;
}

double kyb(double *x, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++)
	if (max < fabs(x[i]))
		max = fabs(x[i]);
	return max;
}