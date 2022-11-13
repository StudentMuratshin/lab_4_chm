#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "fun.h"
using namespace std;

constexpr auto M_PI = 3.14159265358979323846;

double f(double x)
{
	return x - sqrt(x) + 1. / ((3. * x + 1.) * (3. * x + 1.));
}

double df(double x)
{
	return 1 - 1 / (2 * sqrt(x)) - 6 / pow((3 * x + 1), 3);
}

double phi_2(double x)
{
	return pow(x + 1. / ((3. * x + 1.) * (3. * x + 1.)), 2);
}

double phi_1(double x)
{
	return sqrt(x) - 1. / ((3. * x + 1.) * (3. * x + 1.));
}

double phi_11(double x, double y)
{
	return (3 + cos(y) * log(sin(x) + 1)) / 2.;
}

double phi_1_N(double x, double y)
{
	return -sin(y) * cos(x) / (sin(x) + 1);
}

double phi_2_N(double x, double y)
{
	return 2 * y - 3 - cos(y) * log(sin(x) + 1);
}

double phi_1_N_x(pair<double, double> point)
{
	return sin(point.second) * sin(point.first) / (sin(point.first) + 1) + sin(point.second) * cos(point.first) * cos(point.first) / pow((sin(point.first) + 1), 2);
}

double phi_1_N_y(pair<double, double> point)
{
	return -cos(point.second) * cos(point.first) / (sin(point.first) + 1);
}

double phi_2_N_y(pair<double, double> point)
{
	return 2 + sin(point.second) * log(sin(point.first) + 1);
}

double phi_2_N_x(pair<double, double> point)
{
	return -cos(point.second) * cos(point.first) / (sin(point.first) + 1);
}

double Bisection_method(double a, double b, double delta, vector<double>& x)
{
    double c = a;
    while ((b - a) >= delta) {
        c = (a + b) / 2;
        if (abs(f(c)) < std::numeric_limits<double>::epsilon())
		{			
			x.push_back(c);
            return c;
		}

		if (f(a) * f(c) < 0 && f(c) * f(b) < 0) Bisection_method(c, b, delta, x);
        if (f(c) * f(a) < 0)
            b = c;
        else
            a = c;
    }
	x.push_back(c);
    return c;
}

void Chord(double a, double b, double delta, vector<double>& x, string oper)
{
	double c = b, c_prev = b;
	if (oper == "+")
	{
		do
		{
			c_prev = c;
			c = c_prev + f(c_prev) * (c_prev - a) / (f(c_prev) - f(a));
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
	else if (oper == "-")
	{
		do
		{
			c_prev = c;
			c = c_prev - f(c_prev) * (c_prev - a) / (f(c_prev) - f(a));
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
}

void Fixed_point_iteration(double a, double b, double delta, vector<double>& x, int root)
{
	double c = b, c_prev = b;
	if (root == 1)
	{
		do
		{
			c_prev = c;
			c = phi_1(c_prev);
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
	else if (root == 2)
	{
		do
		{
			c_prev = c;
			c = phi_2(c_prev);
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
}

void Newton_method(double a, double b, double delta, std::vector<double>& x, int root)
{
	double c = b, c_prev = b;
	if (root == 1)
	{
		do
		{
			c_prev = c;
			c = c_prev - f(c_prev) / df(c_prev);
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
	else if (root == 2)
	{
		c = a + 1e-10, c_prev = a;
		do
		{
			c_prev = c;
			c = c_prev - f(c_prev) / df(c_prev);
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
}

void Secant_method(double a, double b, double delta, std::vector<double>& x, int root)
{
	double c = a, c_prev = b, c_prev_prev;
	if (root == 1)
	{
		do
		{
			c_prev_prev = c_prev;
			c_prev = c;
			c = c_prev - (c_prev - c_prev_prev) * f(c_prev) / (f(c_prev) - f(c_prev_prev));
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
	else if (root == 2)
	{
		do
		{
			c_prev_prev = c_prev;
			c_prev = c;
			c = c_prev - (c_prev - c_prev_prev) * f(c_prev) / (f(c_prev) - f(c_prev_prev));
		} while (abs(c - c_prev) > delta);
		x.push_back(c);
	}
}

pair<double, double> Minimize(double a, double b, double delta)
{
	double c = a, c_prev = b;
	do
	{
		c_prev = c;
		c = phi_11(M_PI / 2., c_prev);
	} while (abs(c - c_prev) > delta);
	return pair<double, double>(M_PI / 2., c);
}

std::pair<double, double> Minimize_Newton(double a, double b, double delta)
{
	pair <double, double> point = { 1.4 , 1.6 };
	double c[2] = { 1.4 , 1.6 }, c_prev[2];
	double delta_k[2] = { -point.first, -point.second };
	do
	{
		copy(c, c + 2, c_prev);
		c[0] = c_prev[0] + delta_k[0];
		c[1] = c_prev[1] + delta_k[1];
		point = { delta_k[0],delta_k[1] };
		double f[2] = { -point.first, -point.second };
		double M[4] = { phi_1_N_x(point), phi_1_N_y(point), phi_2_N_x(point), phi_2_N_y(point) };
		solve_gauss(2, M, f, delta_k);

	} while (abs(c[0] - c_prev[0]) > delta && abs(c[1] - c_prev[1]) > delta);
	
	return std::pair<double, double>(c[0], c[1]);
}

void solve_gauss(const unsigned int n, double* A, double* b, double* x) {
	// forward
	for (int i = 0; i < n; ++i)
	{
		double pivot = A[i + i * n];
		if (abs(pivot) < 1e-10) {
			// Если пилотный элемент равен нулю
			double max = 0.;
			int max_index = i;
			for (int j = i + 1; j < n; ++j) {
				if (abs(A[i + j * n]) > max) {
					max_index = j;
					max = A[j + j * n];
				}
			}

			for (int k = i; k < n; ++k) {
				std::swap(A[k + max_index * n], A[k + i * n]);
			}
			std::swap(b[i], b[max_index]);
			i--;
			continue;
		}
		for (int j = i; j < n; ++j) {
			A[j + i * n] /= pivot;
		}
		b[i] /= pivot;

		for (int j = i + 1; j < n; ++j) {
			// A[j][i];
			pivot = A[i + j * n];
			for (int k = 0; k < n; ++k) {
				A[k + j * n] -= pivot * A[k + n * i];
			}
			b[j] -= pivot * b[i];
		}
	}
	//backward
	for (int i = n - 1; i >= 0; --i) {
		x[i] = 0;
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[j + i * n] * x[j];
		}
		x[i] += b[i];
	}
}