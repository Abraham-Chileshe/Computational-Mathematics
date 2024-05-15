#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>

using namespace std;

void line(){
    for (int i=0; i<115; i++)
        std::cout<<"-";
    std::cout<<std::endl;
}

double func(double x) {
    return 0.5 * pow(x, 2) + cos(2 * x);
}

double f_derivative(double x, int k) {
    if (k == 1) {
        return x - 2 * sin(2 * x);
    } else if (k == 2) {
        return 1 - 4 * cos(2 * x);
    }
    return pow(-1, (k % 2) + 1) * tgamma(k) / pow(x, k);
}

double middle_rectangular(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += func(a + h * (i + 0.5)) * h;
    }
    return sum;
}

double mr_error(double (*func)(double), double a, double b, int n) {
    double m = 0;
    for (int i = 0; i <= 1000; ++i) {
        m = max(m, abs(f_derivative(a + (b - a) * i / 1000, 2)));
    }
    return m / 24 * pow(b - a, 3) / pow(n, 2);
}

double left_rectangular(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += func(a + h * i) * h;
    }
    return sum;
}

double right_rectangular(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0;
    for (int i = 1; i <= n; ++i) {
        sum += func(a + h * i);
    }
    return sum * h;
}

double trapezoidal(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (func(a) + func(b)) / 2;
    for (int i = 1; i < n; ++i) {
        sum += func(a + h * i);
    }
    return sum * h;
}

double simpson(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0;
    for (int i = 1; i <= n; ++i) {
        sum += func(a + h * (i - 1)) + 4 * func(a + h * (i - 0.5)) + func(a + h * i);
    }
    return sum * h / 6;
}

double l_rect_error(double (*func)(double), double a, double b, int n) {
    double m = 0;
    for (int i = 0; i <= 1000; ++i) {
        m = max(m, abs(f_derivative(a + (b - a) * i / 1000, 1)));
    }
    return m * (b - a) / 2;
}

double r_rect_error(double (*func)(double), double a, double b, int n) {
    double m = 0;
    for (int i = 0; i <= 1000; ++i) {
        m = max(m, abs(f_derivative(a + (b - a) * i / 1000, 1)));
    }
    return m * (b - a) / 2;
}

double trapezoidal_error(double (*func)(double), double a, double b, int n) {
    double m = 0;
    for (int i = 0; i <= 1000; ++i) {
        m = max(m, abs(f_derivative(a + (b - a) * i / 1000, 2)));
    }
    return m / 12 * pow(b - a, 3) / pow(n, 2);
}

double simpson_error(double (*func)(double), double a, double b, int n) {
    double m = 0;
    for (int i = 0; i <= 1000; ++i) {
        m = max(m, abs(f_derivative(a + (b - a) * i / 1000, 4)));
    }
    return m / 2880 * pow(b - a, 5) / pow(n, 4);
}

int main() {
    double a = 0.6, b = 1.1;
    int n = 1;
    double I = 0.124061;

    struct Result {
        vector<int> j;
        vector<int> n;
        vector<double> I_n;
        vector<double> delta_I_n;
        vector<double> relative_I_n;
        vector<double> R_n;
        vector<double> growth;
    } result;

    result.growth.push_back(0);

    for (int i = 0; i < 15; ++i) {
        n *= 2;
        double I_n = middle_rectangular(func, a, b, n);
        result.j.push_back(i + 1);
        result.n.push_back(n);
        result.I_n.push_back(I_n);
        result.delta_I_n.push_back(abs(I - I_n));
        result.relative_I_n.push_back(result.delta_I_n[i] / abs(I) * 100);
        result.R_n.push_back(mr_error(func, a, b, n));
        if (i > 0) {
            result.growth.push_back(result.delta_I_n[i] / result.delta_I_n[i - 1]);
        }
    }
    std::cout<<std::endl<<std::endl;
    cout << "1. Table of values for the central rectangle method:" << endl;
    line();
    cout << setw(10) << "Iteration" << setw(10) << "n" << setw(15) << "I_n" << setw(15) << "delta_I_n" << setw(26) << "Relative Error (%)" << setw(15) << "R_n" << setw(20) << "Growth" << endl;
    line();
    for (size_t i = 0; i < result.j.size(); ++i) {
        cout << setw(10) << result.j[i] << setw(10) << result.n[i] << setw(15) << result.I_n[i] << setw(15) << result.delta_I_n[i] << setw(20) << result.relative_I_n[i] << setw(25) << result.R_n[i] << setw(15) << result.growth[i] << endl;
        line();
    }

    struct Calculate {
        vector<string> method;
        vector<double> I_n;
        vector<double> delta_I_n;
        vector<double> relative_I_n;
        vector<double> R_n;
    } calculate;

    calculate.method = {"Left Rectangles", "Right Rectangles", "Middle Rectangles", "Trapezoids", "Simpson"};

    for (int i = 0; i < 5; ++i) {
        double (*formula)(double (*)(double), double, double, int);
        double (*error)(double (*)(double), double, double, int);
        switch (i) {
            case 0:
                formula = left_rectangular;
                error = l_rect_error;
                break;
            case 1:
                formula = right_rectangular;
                error = r_rect_error;
                break;
            case 2:
                formula = middle_rectangular;
                error = mr_error;
                break;
            case 3:
                formula = trapezoidal;
                error = trapezoidal_error;
                break;
            case 4:
                formula = simpson;
                error = simpson_error;
                break;
        }
        calculate.I_n.push_back(formula(func, a, b, 10000));
        calculate.delta_I_n.push_back(abs(I - calculate.I_n[i]));
        calculate.relative_I_n.push_back(calculate.delta_I_n[i] / abs(I) * 100);
        calculate.R_n.push_back(error(func, a, b, 10000));
    }
    std::cout<<std::endl<<std::endl;
    cout << "2. Table of values for different methods:" << endl;
    line();
    cout << setw(25) << "Method" << setw(15) << "I_n" << setw(15) << "delta_I_n" << setw(20) << "Relative Error (%)" << setw(15) << "R_n" << endl;
    line();
    for (size_t i = 0; i < calculate.method.size(); ++i) {
        cout << setw(25) << calculate.method[i] << setw(15) << calculate.I_n[i] << setw(15) << calculate.delta_I_n[i] << setw(20) << calculate.relative_I_n[i] << setw(15) << calculate.R_n[i] << endl;
    }

    return 0;
}
