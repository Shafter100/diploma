#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <float.h>
#include "omp.h"

#include <shlwapi.h>
#include <windows.h>
#pragma comment (lib, "shell32.lib ")

using namespace std;

int type = 0;

void printm(double*b, int n); // ����� �������
void printm(double**a, int n);
void progonka(double *x, double *a, double *b1, double *c, double *d1, int n); // ��������
double *m3mult(double *a1, double *a2, double *a3, double *b, int n); // ������������ ���������������� "�������" �� ������
double *mmult(double a, double *b, int n); // ������������ ��������� �� ������
double *sum(double *a, double *b, int n); // ����� ��������
void equate(double *a, double *b, int n); // ������������� ��������
double q1(double t); // ����� �����
double q2(double t); // ����� ������
double f1(double x, double t); // ������� ��� �����
double ts(double t, double alpha); // ����������� �����
void zeidel(double**a, double **aij_ii, double *b, double *x, int n); // ����� �������. aij_ii[i][j] = a[i][j] / a[i][i]
void mmult_nloc(double *b1, double **a, double *b, int n); // ������������ ������� ���� ������� ������������� (� ��������� ������ ����������) �� ������

void make_P(double *P, double tau, double t_relax, double **K, double *T_new, double *T_old, int n); // ������������ ������� ��� ������� P
void make_rhs(double *d, double tau, double t_relax, double *f, double *P, double **rhs_dop, double *T, int n); // ������ ����� ��� �����


double *mmult(double **a, double *b, int n);
double **mmult(double **a, double **b, int n);
double kyb(double *x, int n); // ����� ������� ����������

void tauF_cT(double *d, double tau, double *f, double *c1, double *c2, double *c3, double *T, int n); // tau*f + cT (��������� � d)

//������ ������������� (� ���������)
int k = 3;

int main()
{
	setlocale(LC_ALL, "Russian");
	
	double L = 1.; double ro = 1.; double la = 1.; double c = 1.; // ��������� �������
	int n = 100; // ���-�� �������� ���������
	double h = L / n; // ��� �� �������
	n++; // ���������� ����� �� ������� ������ ���������� ���������	
	double tau = 0.01; // ��� �� �������
	double t1 = 3.; // �����
	double m = round(t1 / tau); // ���������� ����� �� �������
	double T0 = 400; // ��������� �������

	double t_relax = 1;

	double p1 = 0.5;
	double p2 = 1 - p1;  // ��������� �������������
	double a = (double)k;

	double *c1 = new double[n];
	double *c2 = new double[n];		
	double *c3 = new double[n];
								// ������� ������� ��� ���������� ������ � � �
	double *k1 = new double[n];
	double *k2 = new double[n];	
	double *k3 = new double[n];

	c1[0] = 0; c2[0] = ro*c*h / 3.; c3[0] = ro*c*h / 6.;
	k1[0] = 0; k2[0] = p1*la / h; k3[0] = -p1*la / h;

	for (int i = 1; i < n - 1; i++)
	{
		c1[i] = ro*c*h / 6.;
		c2[i] = 2. * ro*c*h / 3.;
		c3[i] = ro*c*h / 6.;					// ��������� "�������" � � �
		k1[i] = -p1*la / h;
		k2[i] = 2.*p1*la / h;
		k3[i] = -p1*la / h;
	}

	c1[n - 1] = ro*c*h / 6.; c2[n - 1] = ro*c*h / 3.; c3[n - 1] = 0;
	k1[n - 1] = -p1*la / h; k2[n - 1] = p1*la / h; k3[n - 1] = 0;

	// ��������� Mathematica, ����� ��������� ������� �������������
	cout << "n = " << n << ", a = " << a << endl;
	ShellExecuteA(NULL, "open", "\"A:\\�����\\������\\8 �������\\C++\\�������.nb\"", NULL, NULL, SW_SHOW);
	system("pause");

	//������� �������������

	double **matrixA = new double*[n];
	for (int i = 0; i < n; i++)
		matrixA[i] = new double[n];
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		matrixA[i][j] = 0;

	ifstream input("A:/�����/������/8 �������/C++/matrix_nloc.dat");
	if (input.is_open())
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	{
		input >> matrixA[i][j];
		matrixA[i][j] *= la*p2 / (2 * a*h);
	}	
	else
		cout << "������ ������ �����" << endl;

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

	double **matrixK = new double*[n]; // �������� ������� K ��� ���������� ������� P
	for (int i = 0; i < n; i++)
		matrixK[i] = new double[n];

	double dop = - t_relax*(1 - exp(-tau / t_relax)); // ����������� ��� ������� K+K_nloc ��� ������ ������ ����� ��� ������� �������

	double **rhs_dop= new double*[n]; // �������� C-(K+K_nloc)*dop (����� ����� � dop)
	for (int i = 0; i < n; i++)
		rhs_dop[i] = new double[n];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	{
		rhs_dop[i][j] = matrixA[i][j] * dop;
		matrixK[i][j] = matrixA[i][j];		// ��������� K
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

	// �������� ������� � ��� �����
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

	double *T = new double[n]; // ������ �������
	for (int i = 0; i < n; i++)
		T[i] = T0;

	cout << "\t\t �������� ��� ������: " << endl;
	cout << "1. ���������� ������ �� �������� " << endl << "2. ��������� ������� 2-�� ���� " << endl << "3. ��������� ������� 2-�� � 3-�� ����" << endl << "��� ������ ������� 0" << endl;
	cin >> type;
	if (type != 0) // ���� ��������� 
	{
		if (type == 1) // ���������� ������ �� ��������
		{
			double q1 = 100.; double q2 = 0; // ��������� �������
			double *f = new double[n]; // ������ ������ �����
			f[0] = q1; f[n - 1] = q2;
			for (int i = 1; i < n - 1; i++)			// ��������� ������ ������ �����
				f[i] = 0;

			double *P = new double[n]; // ��������������� ������, ������� ��������� �� ������������ �������
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			//double t1, t2;
			//t1 = omp_get_wtime();

			ofstream fout;
			fout.open("A:/�����/������/8 �������/C++/vecT10.dat");

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				for (int i = 0; i < n; i++)
					fout << T[i] << " ";			// ��� ������ �� ��������� �����
				fout << endl;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);
			}
			//t2 = omp_get_wtime();

			/*ofstream fout;
			fout.open("A:/�����/������/8 �������/C++/vecT1.dat");*/
			//fout.open("F:/������/C++/vecT1.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			//cout << "����� ����������: " << t2 - t1 << endl;

			delete[] f;
			delete[] rhs;
			delete[] T_old;
			delete[] P;
		}

		if (type == 2) // ��������� ������� 2-�� ����
		{
			double *f = new double[n]; // ������ ������ �����
			double t;

			double *P = new double[n]; // ��������������� ������, ������� ��������� �� ������������ �������
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				t = tau*j;
				f[0] = f1(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = f1(L - h / 2, t)*h / 2 - q2(t);
				for (int i = 1; i < n - 1; i++)
					f[i] = (f1((i - 1.0 / 2)*h, t) + f1((i + 1.0 / 2)*h, t))*h / 2;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);
			}

			ofstream fout;
			fout.open("A:/�����/������/8 �������/C++/vecT2.dat");
			//fout.open("F:/������/C++/vecT2.dat");
			for (int i = 0; i < n; i++)
				fout << T[i] << " ";
			fout.close();

			delete[] f;
			delete[] rhs;
			delete[] T_old;
			delete[] P;			
		}

		if (type == 3) // ��������� ������� 2-�� � 3-�� ����
		{
			double *f = new double[n]; // ������ ������ �����
			double t;
			double alpha = 10.; // ����������� ������������� �����������

			matrixA[n - 1][n - 1] += tau*alpha;
			for (int i = 0; i < n; i++)
				aij_ii[i][n - 1] = matrixA[i][n - 1] / matrixA[i][i];
			for (int j = 0; j < n; j++)
				aij_ii[n - 1][j] = matrixA[n - 1][j] / matrixA[n - 1][n - 1];

			double *P = new double[n]; // ��������������� ������, ������� ��������� �� ������������ �������
			for (int i = 0; i < n; i++)
				P[i] = 0;

			double *rhs = new double[n];

			double *T_old = new double[n];

			for (int j = 0; j <= m; j++) // ���� �� �������
			{
				t = tau*j;
				f[0] = f1(h / 2, t)*h / 2 + q1(t);
				f[n - 1] = f1(L - h / 2, t)*h / 2 + alpha*ts(t, alpha);
				for (int i = 1; i < n - 1; i++)
					f[i] = (f1((i - 1.0 / 2)*h, t) + f1((i + 1.0 / 2)*h, t))*h / 2;

				equate(T_old, T, n);

				make_rhs(rhs, tau, t_relax, f, P, rhs_dop, T, n);

				zeidel(matrixA, aij_ii, rhs, T, n);

				make_P(P, tau, t_relax, matrixK, T, T_old, n);
			}
			
			ofstream fout;
			fout.open("A:/�����/������/8 �������/C++/vecT3.dat");
			//fout.open("F:/������/C++/vecT3.dat");
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

//void zeidel(double**a, double **aij_ii, double *b, double *x, int n) // ������������������ ������
//{
//	double  norma_xr;
//
//	double *y, *x0, *xr;
//
//	double epsilon = 1e-3;
//
//	y = new double[n];
//	x0 = new double[n];
//	xr = new double[n];
//
//		for (int i = 0; i < n; i++)
//		{
//			y[i] = b[i] / a[i][i];
//			x0[i] = x[i];
//		}
//
//		int iter = 0;
//		do
//		{
//			for (int i = 0; i < n; i++)
//			{
//				x[i] = y[i];
//
//				for (int j = i + 1; j < n; j++)
//					x[i] -= aij_ii[i][j] * x0[j];
//
//				for (int j = 0; j < i; j++)
//					x[i] -= aij_ii[i][j] * x[j];
//
//				xr[i] = fabs(x0[i] - x[i]);
//				x0[i] = x[i];
//			}
//			iter++;
//
//
//			norma_xr = kyb(xr, n);
//
//		//} while (iter < 20);
//		} while (norma_xr > epsilon);
//
//		/*cout << "\t �������: " << endl;
//		printm(x, n);*/
//
//		/*cout << "������ �� " << iter << " �������� " << endl;
//		system("pause");*/
//
//	delete[] y;
//	delete[] x0;
//	delete[] xr;
//}

void zeidel(double**a, double **aij_ii, double *b, double *x, int n)
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

	/*cout << "������ �� " << iter << " �������� " << endl;
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

void progonka(double *x, double *a, double *b1, double *c, double *d1, int n) // ��������
{
	double *alfa, *beta, *b, *d;
	alfa = new double[n];
	beta = new double[n];
	b = new double[n];		// ��� ����� �����
	d = new double[n];

	for (int i = 0; i < n; i++)
	{
		b[i] = -b1[i];
		d[i] = -d1[i];
	}

	alfa[0] = beta[0] = 0;
	alfa[1] = c[0] / b[0];
	beta[1] = d[0] / b[0];

	for (int i = 1; i <= n - 2; i++)
	{
		alfa[i + 1] = c[i] / (b[i] - a[i] * alfa[i]);
		beta[i + 1] = (d[i] + a[i] * beta[i]) / (b[i] - a[i] * alfa[i]);
	}

	x[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 1]) / (b[n - 1] - a[n - 1] * alfa[n - 1]);

	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = alfa[i + 1] * x[i + 1] + beta[i + 1];
	}

	delete[]alfa;
	delete[]beta;
}

void tauF_cT(double *d, double tau, double *f, double *c1, double *c2, double *c3, double *T, int n) // tau*f + cT
{
	d[0] = tau*f[0] + c2[0] * T[0] + c3[0] * T[1];
	for (int i = 1; i < n - 1; i++)
		d[i] = tau*f[i] + c1[i] * T[i - 1] + c2[i] * T[i] + c3[i] * T[i + 1];
	d[n - 1] = tau*f[n - 1] + c1[n - 1] * T[n - 2] + c2[n - 1] * T[n - 1];

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
		return 4*t*t*exp(-2.0*t);
		//return t;
		//return pow(t,5);
	if (type == 3)
		return t;
		//return 10;
}

double q2(double t)
{
	return 0;
	//return t / exp(1.0);
}

double f1(double x, double t)
{
	if (type == 2)
		//return exp(-pow(x, 2))*(1 - x + 2 * t*(1 + x)*(1 + 2 * (-2 + x)*x));
		return 0;
	if (type == 3)
		return exp(-x*x)*(2 - x + t*(4 + 2 * x*(-3 + 2 * (-2 + x)*x)));
		//return 0; 
}

double ts(double t, double alpha)
{
	//return 0;
	return (t / exp(1.0))*(alpha - 3) / alpha;
}

double *mmult(double **a, double *b, int n)
{
	double *b1;
	b1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		b1[i] = 0;

		for (int j = 0; j < n; j++)

			b1[i] += a[i][j] * b[j];
	}
	return b1;
	delete[] b1;
}

double **mmult(double **a, double **b, int n)
{
	//������� ������������
	double **P;
	P = new double*[n];
	for (int i = 0; i < n; i++)
		P[i] = new double[n];
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		P[i][j] = 0;

	//������������
	for (int k = 0; k < n; k++)
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		P[i][j] += (a[i][k] * b[k][j]);
	return P;

	for (int i = 0; i < n; i++)
		delete[] P[i];
	delete[] P;
}//mmult

double kyb(double *x, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++)
	if (max < fabs(x[i]))
		max = fabs(x[i]);
	return max;
}