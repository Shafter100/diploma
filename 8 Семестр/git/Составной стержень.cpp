#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <float.h>

using namespace std;

typedef double myt;
const myt eps = 1e-7;
int type = 0;

void printm(myt*b, int n); // Вывод вектора
myt *progonka(myt *a, myt *b1, myt *c, myt *d1, int n); // Прогонка
myt *mod_progonka(myt *a, myt *b1, myt *c, myt *d1, int n1, /* Тут разделение */ myt *ap, myt *b1p, myt *cp, myt *d1p, /* Коэффициент */ int n2, myt k);
// Модифицированная прогонка для составного стержня
myt *m3mult(myt *a1, myt *a2, myt *a3, myt *b, int n); // Произведение трехдиагональной "матрицы" на вектор
myt *mmult(myt a, myt *b, int n); // Произведение константы на вектор
myt *sum(myt *a, myt *b, int n); // Сумма векторов
void equate(myt *a, myt *b, int n); // Приравнивание векторов
myt q1(myt t); // Поток слева
myt q2(myt t); // Поток справа
myt func(myt x, myt t); // Функция для теста
myt ts(myt t, myt alpha); // Температура среды

int main()
{
	setlocale(LC_ALL, "Russian");
	
	// Длина стержня
	myt L = 1.; // Параметры стержня
	//myt ro_l = 5400.; myt c_l = 600.;
	//myt la_l = 41.8; // Коэффициент теплопроводности на левой части стержня
	//myt ro_r = 8940.; myt c_r = 485.;
	//myt la_r = 401; // Коэффициент теплопроводности на правой части стержня
	//myt ro_r = 7800.; myt c_r = 460.;
	//myt la_r = 22.4; // Коэффициент теплопроводности на правой части стержня
	myt ro_l = 1.; myt c_l = 1.;
	myt la_l = 1.; // Коэффициент теплопроводности на левой части стержня
	myt ro_r = 1.; myt c_r = 1.;
	myt la_r = 1.; // Коэффициент теплопроводности на правой части стержня
	int n = 100; // Кол-во конечных элементов
	myt h = L / n; // Шаг по стержню
	n++; // Количество узлов на единицу больше количества элементов	
	myt tau = 0.1; // Шаг по времени
	myt t1 = 2.; // Время
	myt m = round(t1 / tau); // Количество шагов по времени
	myt T0 = 1000; // Начальное условие
	// Точка соединения двух частей стержня
	int n1 = 3 * n / 10;
	int n2 = n - n1 + 1;
	//cout << n1 << " " << n2 << endl;

	myt *c1_l = new myt[n1];
	myt *c2_l = new myt[n1];		
	myt *c3_l = new myt[n1];
								// Создаем вектора трёх диагоналей матриц С и К левой части стержня
	myt *k1_l = new myt[n1];
	myt *k2_l = new myt[n1];	
	myt *k3_l = new myt[n1];


	c1_l[0] = 0; c2_l[0] = ro_l*c_l*h / 3; c3_l[0] = ro_l*c_l*h / 6;
	k1_l[0] = 0; k2_l[0] = la_l / h; k3_l[0] = -la_l / h;

	for (int i = 1; i < n1 - 1; i++)
	{
		c1_l[i] = ro_l*c_l*h / 6;
		c2_l[i] = 2 * ro_l*c_l*h / 3;
		c3_l[i] = ro_l*c_l*h / 6; 				// Заполняем "матрицы" С и К левой части стержня
		k1_l[i] = -la_l / h;
		k2_l[i] = 2 * la_l / h;
		k3_l[i] = -la_l / h;
	}

	c1_l[n1 - 1] = ro_l*c_l*h / 6; c2_l[n1 - 1] = ro_l*c_l*h / 3; c3_l[n1 - 1] = 0;
	k1_l[n1 - 1] = -la_l / h; k2_l[n1 - 1] = la_l / h; k3_l[n1 - 1] = 0;


	myt *c1_r = new myt[n2];
	myt *c2_r = new myt[n2];
	myt *c3_r = new myt[n2];
										// Создаем вектора трёх диагоналей матриц С и К правой части стержня
	myt *k1_r = new myt[n2];
	myt *k2_r = new myt[n2];
	myt *k3_r = new myt[n2];


	c1_r[0] = 0; c2_r[0] = ro_r*c_r*h / 3; c3_r[0] = ro_r*c_r*h / 6;
	k1_r[0] = 0; k2_r[0] = la_r / h; k3_r[0] = -la_r / h;

	for (int i = 1; i < n2 - 1; i++)
	{
		c1_r[i] = ro_r*c_r*h / 6;
		c2_r[i] = 2 * ro_r*c_r*h / 3;
		c3_r[i] = ro_r*c_r*h / 6; 				// Заполняем "матрицы" С и К правой части стержня
		k1_r[i] = -la_r / h;
		k2_r[i] = 2 * la_r / h;
		k3_r[i] = -la_r / h;
	}

	c1_r[n2 - 1] = ro_r*c_r*h / 6; c2_r[n2 - 1] = ro_r*c_r*h / 3; c3_r[n2 - 1] = 0;
	k1_r[n2 - 1] = -la_r / h; k2_r[n2 - 1] = la_r / h; k3_r[n2 - 1] = 0;

	/*cout << "C:" << endl;
	printm(c3_l, n1);
	cout << endl;
	printm(c2_l, n1);
	cout << endl;
	printm(c1_l, n1);
	cout << endl;

	printm(c3_r, n2);
	cout << endl;
	printm(c2_r, n2);
	cout << endl;
	printm(c1_r, n2);
	cout << endl;

	cout << "K:" << endl;	
	printm(k3_l, n1);
	cout << endl;
	printm(k2_l, n1);
	cout << endl;
	printm(k1_l, n1);
	cout << endl;

	printm(k3_r, n2);
	cout << endl;
	printm(k2_r, n2);
	cout << endl;
	printm(k1_r, n2);
	cout << endl;*/
	

	myt *T = new myt[n]; // Вектор решения
	for (int i = 0; i < n; i++)
		T[i] = T0;

	cout << "\t\t Выберите тип задачи: " << endl;
	cout << "1. Постоянные потоки " << endl << "2. Граничные условия 2-го рода " << endl << "3. Граничные условия 2-го и 3-го рода" << endl << "Для выхода нажмите 0" << endl;
	cin >> type;
	if (type != 0) // Тело программы 
	{
		if (type == 1) // Постоянные потоки
		{
			myt q1 = 1000.; myt q2 = 0; // Начальные условия
			myt *f = new myt[n]; // Вектор правой части
			f[0] = q1; f[n - 1] = q2;
			for (int i = 1; i < n - 1; i++)			// Заполняем вектор правой части
				f[i] = 0;			
		
			myt *T1 = new myt[n1]; // Решение по частям стержня для цикла по времени (правой части системы по времени)
			myt *T2 = new myt[n2];
			myt *f1 = new myt[n1]; // Правая часть по частям стержня для цикла по времени (правой части системы по времени)
			myt *f2 = new myt[n2];

			myt *cT1 = new myt[n1]; // Произведение CT для разных частей
			myt *cT2 = new myt[n2];

			myt *a1_l = new myt[n1];
			myt *a2_l = new myt[n1];
			myt *a3_l = new myt[n1];

			myt *a1_r = new myt[n2];
			myt *a2_r = new myt[n2];
			myt *a3_r = new myt[n2];

			myt *tauFcT1 = new myt[n1];
			myt *tauFcT2 = new myt[n2];
			myt k;

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				// Коэффициент для решения системы уравнений в месте соединения
				k = -(c1_l[n1 - 1] * T[n1 - 2] + c2_l[n1 - 1] * T[n1 - 1] + c2_r[0] * T[n1 - 1] + c3_r[0] * T[n1]);

				for (int i = 0; i < n1; i++)
				{
					T1[i] = T[i];
					f1[i] = f[i];
				}

				for (int i = 0; i < n2; i++)
				{
					T2[i] = T[i + n1 - 1];
					f2[i] = f[i + n1 - 1];
				}

				cT1 = m3mult(c1_l, c2_l, c3_l, T1, n1);
				cT2 = m3mult(c1_r, c2_r, c3_r, T2, n2);

				tauFcT1 = sum(mmult(tau, f1, n1), cT1, n1);
				tauFcT2 = sum(mmult(tau, f2, n2), cT2, n2);

				equate(a1_l, sum(c1_l, mmult(tau, k1_l, n1), n1), n1);	// Формируем матрицу А1
				equate(a2_l, sum(c2_l, mmult(tau, k2_l, n1), n1), n1);
				equate(a3_l, sum(c3_l, mmult(tau, k3_l, n1), n1), n1);

				equate(a1_r, sum(c1_r, mmult(tau, k1_r, n2), n2), n2);	// Формируем матрицу А2
				equate(a2_r, sum(c2_r, mmult(tau, k2_r, n2), n2), n2);
				equate(a3_r, sum(c3_r, mmult(tau, k3_r, n2), n2), n2);

				equate(T, mod_progonka(
					a1_l, a2_l, a3_l, tauFcT1, n1,
					a1_r, a2_r, a3_r, tauFcT2, n2,
					k), n);
			}
			ofstream fout;
			//fout.open("G:/Курсач/C++/vecT1.dat"); 
			fout.open("A:/Учеба/Курсач/5,6 Семестр/C++/vecT1.dat");
			//fout.open("F:/Курсач/C++/vecT1.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] a1_l;
			delete[] a2_l;
			delete[] a3_l;
			delete[] a1_r;
			delete[] a2_r;
			delete[] a3_r;
			delete[] cT1;
			delete[] cT2;
			delete[] T1;
			delete[] T2;
			delete[] f1;
			delete[] f2;
		}

		if (type == 2) // Граничные условия 2-го рода
		{
			myt *f = new myt[n]; // Вектор правой части
			myt t;

			myt *T1 = new myt[n1]; // Решение по частям стержня для цикла по времени (правой части системы по времени)
			myt *T2 = new myt[n2];
			myt *f1 = new myt[n1]; // Правая часть по частям стержня для цикла по времени (правой части системы по времени)
			myt *f2 = new myt[n2];

			myt *cT1 = new myt[n1]; // Произведение CT для разных частей
			myt *cT2 = new myt[n2];
				
			myt *a1_l = new myt[n1];
			myt *a2_l = new myt[n1];		// Матрица А1
			myt *a3_l = new myt[n1];
			
			myt *a1_r = new myt[n2];
			myt *a2_r = new myt[n2];		// Матрица А2
			myt *a3_r = new myt[n2];

			equate(a1_l, sum(c1_l, mmult(tau, k1_l, n1), n1), n1);	// Формируем матрицу А1
			equate(a2_l, sum(c2_l, mmult(tau, k2_l, n1), n1), n1);
			equate(a3_l, sum(c3_l, mmult(tau, k3_l, n1), n1), n1);

			equate(a1_r, sum(c1_r, mmult(tau, k1_r, n2), n2), n2);	// Формируем матрицу А2
			equate(a2_r, sum(c2_r, mmult(tau, k2_r, n2), n2), n2);
			equate(a3_r, sum(c3_r, mmult(tau, k3_r, n2), n2), n2);

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 - q2(t);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				for (int i = 0; i < n1; i++)
				{
					T1[i] = T[i];
					f1[i] = f[i];
				}

				for (int i = 0; i < n2; i++)
				{
					T2[i] = T[i + n1 - 1];
					f2[i] = f[i + n1 - 1];
				}

				myt k = -(c1_l[n1 - 1] * T[n1 - 2] + c2_l[n1 - 1] * T[n1 - 1] + c2_r[0] * T[n1 - 1] + c3_r[0] * T[n1] + tau*(f1[n1 - 1]));

				cT1 = m3mult(c1_l, c2_l, c3_l, T1, n1);
				cT2 = m3mult(c1_r, c2_r, c3_r, T2, n2);

				equate(T, mod_progonka(
					a1_l, a2_l, a3_l, sum(mmult(tau, f1, n1), cT1, n1), n1,
					a1_r, a2_r, a3_r, sum(mmult(tau, f2, n2), cT2, n2), n2, 
					k), n);
			}

			ofstream fout;
			fout.open("A:/Учеба/Курсач/5,6 Семестр/C++/vecT2.dat");
			//fout.open("F:/Курсач/C++/vecT2.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] a1_l;
			delete[] a2_l;
			delete[] a3_l;
			delete[] a1_r;
			delete[] a2_r;
			delete[] a3_r;
			delete[] cT1;
			delete[] cT2;
			delete[] T1;
			delete[] T2;
			delete[] f1;
			delete[] f2;
		}
		

		if (type == 3) // Граничные условия 2-го и 3-го рода
		{
			myt *f = new myt[n]; // Вектор правой части
			myt t;

			myt *T1 = new myt[n1]; // Решение по частям стержня для цикла по времени (правой части системы по времени)
			myt *T2 = new myt[n2];
			myt *f1 = new myt[n1]; // Правая часть по частям стержня для цикла по времени (правой части системы по времени)
			myt *f2 = new myt[n2];

			myt *cT1 = new myt[n1]; // Произведение CT для разных частей
			myt *cT2 = new myt[n2];

			myt alpha = 10; // Коэффициент конвективного теплообмена

			myt *a1_l = new myt[n1];
			myt *a2_l = new myt[n1];		// Матрица А1
			myt *a3_l = new myt[n1];

			myt *a1_r = new myt[n2];
			myt *a2_r = new myt[n2];		// Матрица А2
			myt *a3_r = new myt[n2];

			equate(a1_l, sum(c1_l, mmult(tau, k1_l, n1), n1), n1);	// Формируем матрицу А1
			equate(a2_l, sum(c2_l, mmult(tau, k2_l, n1), n1), n1);
			equate(a3_l, sum(c3_l, mmult(tau, k3_l, n1), n1), n1);

			equate(a1_r, sum(c1_r, mmult(tau, k1_r, n2), n2), n2);	// Формируем матрицу А2
			equate(a2_r, sum(c2_r, mmult(tau, k2_r, n2), n2), n2);
			equate(a3_r, sum(c3_r, mmult(tau, k3_r, n2), n2), n2);

			a2_r[n2 - 1] += tau*alpha;

			for (int j = 0; j <= m; j++) // Цикл по времени
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 + alpha*ts(t, alpha);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				for (int i = 0; i < n1; i++)
				{
					T1[i] = T[i];
					f1[i] = f[i];
				}

				for (int i = 0; i < n2; i++)
				{
					T2[i] = T[i + n1 - 1];
					f2[i] = f[i + n1 - 1];
				}

				myt k = -(c1_l[n1 - 1] * T[n1 - 2] + c2_l[n1 - 1] * T[n1 - 1] + c2_r[0] * T[n1 - 1] + c3_r[0] * T[n1] + tau*(f1[n1 - 1]));

				cT1 = m3mult(c1_l, c2_l, c3_l, T1, n1);
				cT2 = m3mult(c1_r, c2_r, c3_r, T2, n2);

				equate(T, mod_progonka(
					a1_l, a2_l, a3_l, sum(mmult(tau, f1, n1), cT1, n1), n1,
					a1_r, a2_r, a3_r, sum(mmult(tau, f2, n2), cT2, n2), n2,
					k), n);
			}

			ofstream fout;
			fout.open("A:/Учеба/Курсач/5,6 Семестр/C++/vecT3.dat");
			//fout.open("F:/Курсач/C++/vecT3.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();
			delete[] a1_l;
			delete[] a2_l;
			delete[] a3_l;
			delete[] a1_r;
			delete[] a2_r;
			delete[] a3_r;
			delete[] cT1;
			delete[] cT2;
			delete[] T1;
			delete[] T2;
			delete[] f1;
			delete[] f2;
		}
		printm(T, n);
	}

	system("pause");
	delete[] c1_l;
	delete[] c2_l;
	delete[] c3_l;
	delete[] k1_l;
	delete[] k2_l;
	delete[] k3_l;
	delete[] c1_r;
	delete[] c2_r;
	delete[] c3_r;
	delete[] k1_r;
	delete[] k2_r;
	delete[] k3_r;
}//main


myt *progonka(myt *a, myt *b1, myt *c, myt *d1, int n)
{
	
	myt *x, *alfa, *beta, *b, *d;
	x = new myt[n];
	alfa = new myt[n];
	beta = new myt[n];
	b = new myt[n];		// Для смены знака
	d = new myt[n];


	for (int i = 0; i < n; i++)
	{
		b[i] = -b1[i];
		d[i] = -d1[i];
	}

	alfa[0] = c[0] / b[0];
	beta[0] = d[0] / b[0];

	for (int i = 0; i <= n - 2; i++)
	{
		alfa[i + 1] = c[i] / (b[i] - a[i] * alfa[i]);
		beta[i + 1] = (d[i] + a[i] * beta[i]) / (b[i] - a[i] * alfa[i]);
	}

	x[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 1]) / (b[n - 1] - a[n - 1] * alfa[n - 1]);

	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = alfa[i + 1] * x[i + 1] + beta[i + 1];
	}

	return x;

	delete[]x;
	delete[]b;
	delete[]d;
	delete[]alfa;
	delete[]beta;
}

myt *mod_progonka(myt *a, myt *b1, myt *c, myt *d1, int n1, /* Тут разделение */ myt *ap, myt *b1p, myt *cp, myt *d1p, int n2, /* Коэффициент */ myt k)
{
	myt *x_res = new myt[n1 + n2 - 1]; // Результирующий вектор

	myt *x, *alfa, *beta, *b, *d;
	x = new myt[n1];
	alfa = new myt[n1];
	beta = new myt[n1];
	b = new myt[n1];		// Для смены знака
	d = new myt[n1];

	myt *xp, *delta, *gamma, *bp, *dp;
	xp = new myt[n2];
	delta = new myt[n2];
	gamma = new myt[n2];
	bp = new myt[n2];		// Для смены знака
	dp = new myt[n2];

	for (int i = 0; i < n1; i++) // Левая
	{
		b[i] = -b1[i];
		d[i] = -d1[i];					
	}

	for (int i = 0; i < n2; i++) // Правая
	{
		bp[i] = -b1p[i];
		dp[i] = -d1p[i];
	}

	
	//printm(ap, n2);
	//cout << endl;
	//printm(bp, n2);
	//cout << endl; 
	//printm(cp, n2);
	//cout << endl;
	//printm(dp, n2);
	//cout << endl;
	

	alfa[0] = beta[0] = 0;
	alfa[1] = c[0] / b[0];
	beta[1] = d[0] / b[0];

	for (int i = 1; i <= n1 - 2; i++) // Левая
	{
		alfa[i + 1] = c[i] / (b[i] - a[i] * alfa[i]);
		beta[i + 1] = (d[i] + a[i] * beta[i]) / (b[i] - a[i] * alfa[i]);			
	}
		
	/*printm(alfa, n1);
	cout << endl;
	printm(beta, n1);
	cout << endl;

	system("pause");*/

	delta[n2 - 1] = gamma[n2 - 1] = 0;
	delta[n2 - 2] = ap[n2 - 1] / bp[n2 - 1];
	gamma[n2 - 2] = dp[n2 - 1] / bp[n2 - 1];

	for (int i = n2 - 2; i >= 1; i--) // Правая
	{
		delta[i - 1] = ap[i] / (bp[i] - cp[i] * delta[i]);
		gamma[i - 1] = (dp[i] + cp[i] * gamma[i]) / (bp[i] - cp[i] * delta[i]);
	}

	
	/*printm(delta, n2);
	cout << endl;
	printm(gamma, n2);
	cout << endl;
	
	system("pause");*/

	// Из решения системы на соединении двух частей 
	//printm(beta, n1);
	//printm(alfa,n1);
	//printm(delta, n2);
	//printm(gamma, n2);

	//cout << k << endl;

	/*d[n1 - 1] = ((b[n1 - 1] - a[n1 - 1] * alfa[n1 - 1])*(k + cp[0] * gamma[0]) - (bp[0] - cp[0] * delta[0])*a[n1 - 1] * beta[n1 - 1]) / (bp[0] - cp[0] * delta[0] - b[n1 - 1] + a[n1 - 1] * alfa[n1 - 1]);
	dp[0] = d[n1 - 1] + k;*/

	d[n1 - 1] = ((b[n1 - 1] - a[n1 - 1] * alfa[n1 - 1])*(k + cp[0] * gamma[0]) - (bp[0] - cp[0] * delta[0])*a[n1 - 1] * beta[n1 - 1]) / (bp[0] + b[n1 - 1] - cp[0] * delta[0] - a[n1 - 1] * alfa[n1 - 1]);

	dp[0] = k - d[n1 - 1];

	x[n1 - 1] = (d[n1 - 1] + a[n1 - 1] * beta[n1 - 1]) / (b[n1 - 1] - a[n1 - 1] * alfa[n1 - 1]);

	xp[0] = (dp[0] + cp[0] * gamma[0]) / (bp[0] - cp[0] * delta[0]);

	for (int i = n1 - 2; i >= 0; i--)
	{
		x[i] = alfa[i + 1] * x[i + 1] + beta[i + 1];
	}

	for (int i = 1; i < n2; i++)
	{
		xp[i] = delta[i - 1] * xp[i - 1] + gamma[i - 1];
	}

	// Собираем ответ
	for (int i = 0; i < n1; i++)
		x_res[i] = x[i];

	for (int i = 0; i < n2; i++)
		x_res[i + n1 - 1] = xp[i];

	return x_res;

	delete[] x;
	delete[] b;
	delete[] d;
	delete[] alfa;
	delete[] beta;

	delete[] xp;
	delete[] bp;
	delete[] dp;
	delete[] delta;
	delete[] gamma;

	delete[] x_res;
}

myt *m3mult(myt *a1, myt *a2, myt *a3, myt *b, int n)
{
	myt *b1;
	b1 = new myt[n];
	b1[0] = a2[0] * b[0] + a3[0] * b[1];
	for (int i = 1; i < n - 1; i++)
		b1[i] = a1[i] * b[i - 1] + a2[i] * b[i] + a3[i] * b[i + 1];
	b1[n - 1] = a1[n - 1] * b[n - 2] + a2[n - 1] * b[n - 1];
	return b1;
	delete[] b1;
}

myt *mmult(myt a, myt *b, int n)
{
	myt *b1;
	b1 = new myt[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = a*b[i];
	}
	return b1;
	delete[] b1;
}

myt *sum(myt *a, myt *b, int n)
{
	myt *b1;
	b1 = new myt[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = a[i] + b[i];
	}
	return b1;
	delete[] b1;
}

void printm(myt*b, int n)
{
	for (int i = 0; i < n; i++)
		cout << b[i] << " ";
	cout << endl;
}

void equate(myt *a, myt *b, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = b[i];
}

myt q1(myt t)
{
	if (type == 2)
		return t;
		//return pow(t, 5);
	if (type == 3)
		return t;
		/*return pow(t, 10);*/
		//return 100;
}

myt q2(myt t)
{
	return t/exp(1.0);
}

myt func(myt x, myt t)
{
	if (type == 2)
		return exp(-x*x)*(1 - x + 2*t*(1 + x)*(1 + 2*x*(-2 + x)));
		//return 0;
	if (type == 3)
		return exp(-x*x)*(2 - x + t*(4 + 2 * x*(-3 + 2 * (-2 + x)*x)));
		//return 0; 
}

myt ts(myt t, myt alpha)
{
	//return 0;
	return (t / exp(1.0))*(alpha - 3) / alpha;
}