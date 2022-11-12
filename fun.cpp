#include <cmath>
#include <string>
#include <vector>
#include "fun.h"
using namespace std;



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
