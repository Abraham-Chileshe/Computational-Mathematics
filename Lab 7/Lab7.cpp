#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

double function(double x) {
    return log(0.5 * x) - 0.5 * cos(x);
}

double derivative_function(double x) {
    return 1 / x + 0.5 * sin(x);
}

void print_table(const std::vector<std::pair<int, double>>& table) {
    std::cout << std::setw(10) << "Итерация" << std::setw(20) << "Значение x" << std::endl;
    for (const auto& row : table) {
        std::cout << std::setw(10) << row.first << std::setw(20) << row.second << std::endl;
    }
}

double chord_method(double eps, std::vector<double> lst) {
    std::vector<std::pair<int, double>> table;
    int counter = 1;
    double b = lst[1];
    double x_prev = lst[0];
    double x_next = x_prev - (b - x_prev) * function(x_prev) / (function(b) - function(x_prev));
    table.push_back({counter, x_next});

    while (fabs(x_next - x_prev) > eps) {
        x_prev = x_next;
        x_next = x_prev - (b - x_prev) * function(x_prev) / (function(b) - function(x_prev));
        counter++;
        std::cout << x_next << std::endl;
        table.push_back({counter, x_next});
    }

    print_table(table);
    return x_next;
}

double newton_method(double eps, double x0) {
    std::vector<std::pair<int, double>> table;
    int count = 1;
    double x_prev = x0;
    double x_next = x_prev - function(x_prev) / derivative_function(x_prev);
    table.push_back({count, x_next});

    while (fabs(x_next - x_prev) > eps) {
        count++;
        x_prev = x_next;
        x_next = x_prev - function(x_prev) / derivative_function(x_prev);
        std::cout << x_next << std::endl;
        table.push_back({count, x_next});
    }

    print_table(table);
    return x_next;
}

double bisection_method(double eps, std::vector<double> lst) {
    double l = lst[0];
    double r = lst[1];
    double c = 0;
    int count = 1;

    double x_prev = r;
    double x_next = c;

    while (fabs(x_next - x_prev) > eps) {
        c = (r + l) / 2;

        if (function(c) * function(l) > 0) {
            l = c;
        } else if (function(c) * function(r) > 0) {
            r = c;
        }

        x_prev = x_next;
        x_next = c;
        std::cout << count << " " << x_next << std::endl;
        count++;
    }

    return x_next;
}

int main() {
    std::cout << "Метод хорд:" << std::endl;
    chord_method(0.0000001, {0.5, 2.0});

    std::cout << "\nМетод Ньютона:" << std::endl;
    newton_method(0.0000001, 2.0);

    std::cout << "\nМетод бисекции:" << std::endl;
    bisection_method(0.0000001, {0.5, 2.0});

    return 0;
}
