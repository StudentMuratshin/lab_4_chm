#include <iostream>
#include <vector>
#include "fun.h"

using namespace std;

constexpr auto delta = 1e-6;

int main()
{
    double a = 0., b = 1.;
    vector<double> x;

    cout << "______________Bisection method______________" << endl;
    Bisection_method(a, b, delta, x);
    for (auto t : x) cout << t << " " << endl;
    x.clear();

    cout << "______________Chord method______________" << endl;
    Chord(a, b, delta, x, "+");
    a = 0., b = 0.5;
    Chord(a, b, delta, x, "-");
    for (auto t : x) cout << t << " " << endl;
    x.clear();

    cout << "______________Fixed-point iteration______________" << endl;
    a = 0., b = 1;
    Fixed_point_iteration(a, b, delta, x, 1);
    a = 0., b = 0.5;
    Fixed_point_iteration(a, b, delta, x, 2);
    for (auto t : x) cout << t << " " << endl;
    x.clear();

    cout << "______________Newton method______________" << endl;
    a = 0., b = 1;
    Newton_method(a, b, delta, x, 1);
    a = 0., b = 0.5;
    Newton_method(a, b, delta, x, 2);
    for (auto t : x) cout << t << " " << endl;
    x.clear();

    cout << "______________Secant method______________" << endl;
    a = 0., b = 1;
    Secant_method(a, b, delta, x, 1);
    a = 0., b = 0.5;
    Secant_method(a, b, delta, x, 2);
    for (auto t : x) cout << t << " " << endl;
    x.clear();
}