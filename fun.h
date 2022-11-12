#pragma once
#include <string>
#include <vector>

double f(double x);
double df(double x);
double phi_1(double x);
double phi_2(double x);
double Bisection_method(double a, double b, double delta, std::vector<double> &x);
void Chord(double a, double b, double delta, std::vector<double>& x, std:: string oper);
void Fixed_point_iteration(double a, double b, double delta, std::vector<double>& x, int root);
void Newton_method(double a, double b, double delta, std::vector<double>& x, int root);
void Secant_method(double a, double b, double delta, std::vector<double>& x, int root);