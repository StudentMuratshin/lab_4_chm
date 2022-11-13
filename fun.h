#pragma once
#include <string>
#include <vector>

double f(double x);
double df(double x);
double phi_1(double x);
double phi_2(double x);
double phi_11(double x, double y);
double phi_1_N(double x, double y);
double phi_2_N(double x, double y);

double phi_1_N_x(std::pair<double, double> point);
double phi_2_N_x(std::pair<double, double> point);
double phi_1_N_y(std::pair<double, double> point);
double phi_2_N_y(std::pair<double, double> point);

double Bisection_method(double a, double b, double delta, std::vector<double> &x);
void Chord(double a, double b, double delta, std::vector<double>& x, std:: string oper);
void Fixed_point_iteration(double a, double b, double delta, std::vector<double>& x, int root);
void Newton_method(double a, double b, double delta, std::vector<double>& x, int root);
void Secant_method(double a, double b, double delta, std::vector<double>& x, int root);
std::pair<double, double> Minimize(double a, double b,double delta);
std::pair<double, double> Minimize_Newton(double a, double b, double delta);
void solve_gauss(const unsigned int n, double* A, double* b, double* x);