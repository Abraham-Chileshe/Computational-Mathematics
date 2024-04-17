#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

double newton_minus_param(double t, int n) {
    double a = 1;
    for (int i = 0; i < n; i++) {
        a *= (t - i);
    }
    a /= tgamma(n + 1);
    return a;
}

double newton_plus_param(double t, int n) {
    double a = 1;
    for (int i = 0; i < n; i++) {
        a *= (t + i);
    }
    a /= tgamma(n + 1);
    return a;
}

double gauss1_minus_param(double t, int n) {
    double a = 1;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 1 || i == 0) {
            a *= (t - i);
        } else {
            a *= (t + i - 1);
        }
    }
    a /= tgamma(n + 1);
    return a;
}

double gauss2_plus_param(double t, int n) {
    double a = 1;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 1 || i == 0) {
            a *= (t + i);
        } else {
            a *= (t - i + 1);
        }
    }
    a /= tgamma(n + 1);
    return a;
}

double insert_gauss1_polynomial(double t, int n, const vector<vector<double>>& mass) {
    double Px = 0;
    int j = 5;
    for (int i = 0; i < n; i++) {
        Px += mass[i][j] * gauss1_minus_param(t, i);
        if (i % 2 != 0) {
            j--;
        }
    }
    return Px;
}

double insert_gauss2_polynomial(double t, int n, const vector<vector<double>>& mass) {
    double Px2 = 0;
    int j = 5;
    for (int i = 0; i < n; i++) {
        Px2 += mass[i][j] * gauss2_plus_param(t, i);
        if (i % 2 == 0) {
            j--;
        }
    }
    return Px2;
}

double insert_newton1_polynomial(double t, int n, const vector<vector<double>>& mass) {
    double Px = 0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        Px += mass[i][j] * newton_minus_param(t, i);
    }
    return Px;
}

double insert_newton2_polynomial(double t, int n, const vector<vector<double>>& mass) {
    double Px2 = 0;
    for (int i = 0; i < n; i++) {
        int j = n - i - 1;
        Px2 += mass[i][j] * newton_plus_param(t, i);
    }
    return Px2;
}

int main() {
    const double a = 0.6;
    const double b = 1.1;
    const int n = 11;
    const double h = (b - a) / 10;

    vector<double> x_list(n);
    vector<double> y_list(n);

    for (int i = 0; i < n; i++) {
        double xi = a + i * h;
        x_list[i] = xi;
        y_list[i] = 0.5 * pow(xi, 2) - cos(2 * xi);
    }

    cout << setw(5) << "No." << setw(18) << "x" << setw(18) << "y(x)" << endl;
    cout<<"----------------------------------------------"<<endl;
    for (int i = 0; i < n; i++) {
        cout << setw(5) << i << setw(18) << x_list[i] << setw(18) << y_list[i] << endl;
    }
    cout << endl<<endl;

    vector<vector<double>> list_diffs = {y_list};

    while (list_diffs.back().size() != 1) {
        vector<double> lis;
        for (int i = 0; i < static_cast<int>(list_diffs.back().size()) - 1; i++) {
            lis.push_back(list_diffs.back()[i + 1] - list_diffs.back()[i]);
        }
        list_diffs.push_back(lis);
    }

    vector<vector<double>> list_to_table = list_diffs;
    int max_length = 0;
    for (const auto& lst : list_to_table) {
        max_length = max(max_length, static_cast<int>(lst.size()));
    }

    for (auto& lst : list_to_table) {
        while (lst.size() < static_cast<size_t>(max_length)) {
            lst.push_back(0); // Use 0 or NaN instead of ""
        }
    }

    cout << setw(8) << "No." << setw(18) << "Value 1" << setw(18) << "Value 2" << setw(18) << "Value 3"
     << setw(18) << "Value 4" << setw(18) << "Value 5" << setw(18) << "Value 6" << endl;
    cout<<"-------------------------------------------------------------------------------------------------------------------------"<<endl;
    for (int i = 0; i < static_cast<int>(list_to_table.size()); i++) {
        cout << setw(8) << i;
        for (int j = 0; j < 6; j++) {
            cout << setw(18) << list_to_table[i][j];
        }
        cout << endl<<endl;
    }
    cout << endl<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout << setw(8) << "No."<< setw(18) << "Value 7" << setw(18) << "Value 8" << setw(18) << "Value 9"
         << setw(18) << "Value 10" << setw(18) << "Value 11" << endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------"<<endl;

    for (int i = 0; i < static_cast<int>(list_to_table.size()); i++) {
        cout << setw(8) << i;
        for (int j = 6; j < static_cast<int>(list_to_table[i].size()); j++) {
            cout << setw(18) << list_to_table[i][j];
        }
        cout << endl;
    }


    const double x_star2 = 0.62;
    const double x_star3 = 1.07;
    const double x_star4 = 0.83;

    double t = min(abs(x_list[0] - x_star2), abs(x_list[1] - x_star2)) / h;

    cout<<endl;
    cout << "Newton 1: " << insert_newton1_polynomial(t, n, list_diffs) << endl;
    cout << "R_N1: " << insert_newton1_polynomial(t, n, list_diffs) - (0.5 * pow(x_star2, 2) - cos(2 * x_star2)) << endl;

    t = -1 * (x_list[n - 1] - x_star3) / h;
    cout << "Newton 2: " << insert_newton2_polynomial(t, n, list_diffs) << endl;
    cout << "R_N2: " << insert_newton2_polynomial(t, n, list_diffs) - (0.5 * pow(x_star3, 2) - cos(2 * x_star3)) << endl;

    t = min(abs(x_list[0] - x_star4), abs(x_list[1] - x_star4)) / h;
    cout << "Gauss 1: " << insert_gauss1_polynomial(t, n, list_diffs) << endl;
    cout << "R_G1: " << insert_gauss1_polynomial(t, n, list_diffs) - (0.5 * pow(x_star4, 2) - cos(2 * x_star4)) << endl;

    t = -1 * (x_list[n - 1] - x_star4) / h;
    cout << "Gauss 2: " << insert_gauss2_polynomial(t, n, list_diffs) << endl;
    cout << "R_G2: " << insert_gauss2_polynomial(t, n, list_diffs) - (0.5 * pow(x_star4, 2) - cos(2 * x_star4)) << endl;

    double R_n = insert_newton1_polynomial(t, n, list_diffs);
    double min_Rn = min(R_n, insert_newton1_polynomial(t, n, list_diffs));
    double max_Rn = max(R_n, insert_newton1_polynomial(t, n, list_diffs));
    cout << "Minimum Rn on interval: " << min_Rn << endl;
    cout << "Maximum Rn on interval: " << max_Rn << endl;

    return 0;
}
