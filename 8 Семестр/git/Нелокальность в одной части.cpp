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

void printm(double*b, int n); // ����� �������
void printm(double**a, int n);
double q1(double t); // ����� �����
double q2(double t); // ����� ������
double func(double x, double t); // ������� ��� �����
double ts(double t, double alpha); // ����������� �����
void tauF_cT(double *d, double tau, double *f, double *c1, double *c2, double *c3, double *T, int n); // tau*f + cT

void zeidel(double**a, double **aij_ii, double*b, double*x, int n); // ����� �������. aij_ii[i][j] = a[i][j] / a[i][i]
double kyb(double *x, int n); // ����� ������� ����������

//������ ������������� (� ���������)
int k = 3;

int main()
{
	setlocale(LC_ALL, "Russian");
	
	double L = 1.; // ��������� �������

	double ro_l = 5400.; double c_l = 600.; // ������ ������
	double la_l = 41.8; // ����������� ���������������� �� ����� ����� �������

	//double ro_r = 8940.; double c_r = 485.; // ����
	//double la_r = 401; // ����������� ���������������� �� ������ ����� �������

	double ro_r = 7800.; double c_r = 460.; // �����
	double la_r = 22.4; // ����������� ���������������� �� ������ ����� �������
		
	//double ro_l = 1.; double c_l = 1.;
	//double la_l = 1.; // ����������� ���������������� �� ����� ����� �������
	//double ro_r = 1.; double c_r = 1.;
	//double la_r = 1.; // ����������� ���������������� �� ������ ����� �������
	int n = 100; // ���-�� �������� ���������
	double h = L / n; // ��� �� �������
	n++; // ���������� ����� �� ������� ������ ���������� ���������	
	double tau = 0.01; // ��� �� �������
	double t1 = 1e3; // �����
	double m = round(t1 / tau); // ���������� ����� �� �������
	double T0 = 1000; // ��������� �������
	int n1 = 5 * n / 10; // ����� ���������� ���� ������ �������
	int n2 = n - n1 + 1;
	//cout << n1 << " " << n2 << endl;

	double p1 = 0.75;
	double p2 = 1 - p1;  // ��������� ������������� ��� ����� ����� �������
	double a = (double) k;

	double *c1 = new double[n];
	double *c2 = new double[n];		
	double *c3 = new double[n];
								// ������� ������� ��� ���������� ������ � � � 
	double *k1 = new double[n];
	double *k2 = new double[n];	
	double *k3 = new double[n];

	// ������� ������������ ������ ��� ����� ����� �������

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

	// � ����� ����������
	c1[n1 - 1] = ro_l*c_l*h / 6; c2[n1 - 1] = ro_l*c_l*h / 3 + ro_r*c_r*h / 3; c3[n1 - 1] = ro_r*c_r*h / 6;;
	k1[n1 - 1] = -p1*la_l / h; k2[n1 - 1] = p1*la_l / h + la_r / h; k3[n1 - 1] = -la_r / h;

	// ������� ������������ ��� ������ ����� �������
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

	// ��������� Mathematica, ����� ��������� ������� �������������
	cout << "n1 = " << n1 << ", a = " << a << endl;
	//ShellExecuteA(NULL, "open", "\"A:\\�����\\������\\7 �������\\C++\\�������.nb\"", NULL, NULL, SW_SHOW); // ����� ��������� ����, ��� ����� �������.nb
	system("pause");

	// ������� �������������, ������� ������������ ������ �������� A ��� �������

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

	ifstream input("A:/�����/������/7 �������/C++/matrix_nloc.dat");
	if (input.is_open())
	for (int i = 0; i < n1; i++)
	for (int j = 0; j < n1; j++)
	{
		input >> matrixA[i][j];
		matrixA[i][j] *= la_l*p2 / (2 * a*h);
	}
	else
		cout << "������ ������ �����" << endl;

	//printm(matrixA, n1);
	input.close();

	// ��������� ������� K
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
		
	// �������� ������� � ��� �����
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		matrixA[i][j] *= tau; // tau*K

	// tau*K + C
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


	/*cout << "C:" << endl;
	printm(c3, n);
	cout << endl;
	printm(c2, n);
	cout << endl;
	printm(c1, n);
	cout << endl;

	cout << "K:" << endl;	
	printm(k3, n);
	cout << endl;
	printm(k2, n);
	cout << endl;
	printm(k1, n);
	cout << endl;*/


	double *T = new double[n]; // ������ �������
	for (int i = 0; i < n; i++)
		T[i] = T0;

	cout << "\t\t �������� ��� ������: " << endl;
	cout << "1. ���������� ������ " << endl << "2. ��������� ������� 2-�� ���� " << endl << "3. ��������� ������� 2-�� � 3-�� ����" << endl << "��� ������ ������� 0" << endl;
	cin >> type;
	if (type != 0) // ���� ��������� 
	{
		if (type == 1) // ���������� ������
		{
			double q1 = 1e7; double q2 = 0; // ��������� �������
			double *f = new double[n]; // ������ ������ �����
			f[0] = q1; f[n - 1] = q2;
			for (int i = 1; i < n - 1; i++)			// ��������� ������ ������ �����
				f[i] = 0;			

			double *tauFcT = new double[n];

			double t1, t2;

			t1 = omp_get_wtime();

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				tauF_cT(tauFcT, tau, f, c1, c2, c3, T, n);

				zeidel(matrixA, aij_ii, tauFcT, T, n);
			}
			t2 = omp_get_wtime();
			cout << "�����: " << t2 - t1 << " c" << endl;

			ofstream fout;
			//fout.open("G:/������/C++/vecT1.dat"); 
			fout.open("A:/�����/������/7 �������/C++/vecT1.dat");
			//fout.open("F:/������/C++/vecT1.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] f;
			delete[] tauFcT;
		}

		if (type == 2) // ��������� ������� 2-�� ����
		{
			double *f = new double[n]; // ������ ������ �����
			double t;

			double *tauFcT = new double[n];

			double t1, t2;

			t1 = omp_get_wtime();

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 - q2(t);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				tauF_cT(tauFcT, tau, f, c1, c2, c3, T, n);

				zeidel(matrixA, aij_ii, tauFcT, T, n);
			}

			t2 = omp_get_wtime();
			cout << "�����: " << t2 - t1 << " c" << endl;


			ofstream fout;
			fout.open("A:/�����/������/7 �������/C++/vecT2.dat");
			//fout.open("F:/������/C++/vecT2.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();
			
			delete[] f;
			delete[] tauFcT;
		}
		

		if (type == 3) // ��������� ������� 2-�� � 3-�� ����
		{
			double *f = new double[n]; // ������ ������ �����
			double t;

			double alpha = 10; // ����������� ������������� �����������

			matrixA[n - 1][n - 1] += tau*alpha;

			double *tauFcT = new double[n];

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				t = tau*j;
				f[0] = func(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = func(L - h / 2, t)*h / 2 + alpha*ts(t, alpha);
				for (int i = 1; i < n - 1; i++)
					f[i] = (func((i - 1.0 / 2)*h, t) + func((i + 1.0 / 2)*h, t))*h / 2;

				tauF_cT(tauFcT, tau, f, c1, c2, c3, T, n);

				zeidel(matrixA, aij_ii, tauFcT, T, n);
			}

			ofstream fout;
			fout.open("A:/�����/������/7 �������/C++/vecT3.dat");
			//fout.open("F:/������/C++/vecT3.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] f;
			delete[] tauFcT;
		}
		//printm(T, n);
	}

	system("pause");
	for (int i = 0; i < n; i++)
	{
		delete[] matrixA[i];
		delete[] aij_ii[i];
	}
	delete[] matrixA;
	delete[] aij_ii;

	delete[] c1;
	delete[] c2;
	delete[] c3;
	delete[] k1;
	delete[] k2;
	delete[] k3;
}//main

void tauF_cT(double *d, double tau, double *f, double *c1, double *c2, double *c3, double *T, int n) // tau*f + cT
{
	d[0] = tau*f[0] + c2[0] * T[0] + c3[0] * T[1];
	for (int i = 1; i < n - 1; i++)
		d[i] = tau*f[i] + c1[i] * T[i - 1] + c2[i] * T[i] + c3[i] * T[i + 1];
	d[n - 1] = tau*f[n - 1] + c1[n - 1] * T[n - 2] + c2[n - 1] * T[n - 1];

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

//������ �����
double func(double x, double t)
{
	if (type == 2)
		return exp(-x*x)*(1 - x + 2*t*(1 + x)*(1 + 2*x*(-2 + x)));
	//return 0;
	if (type == 3)
		//return exp(-x*x)*(2 - x + t*(4 + 2 * x*(-3 + 2 * (-2 + x)*x)));
		return 0; 
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

	//double r0 = kyb(mmult_min(a, x0, b, n), n);

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

		//norma_xr = kyb(mmult_min(a, x, b, n), n) / r0;
		norma_xr = kyb(xr, n);
	//} while (iter < 100);
	} while (norma_xr > epsilon);
	//} while (norma_xr > (1 - norma_c) / norma_c * epsilon);

	/*cout << "������ �� " << iter << " �������� " << endl;
	system("pause");*/

	delete[] y;
	delete[] x0;
	delete[] xr;
}