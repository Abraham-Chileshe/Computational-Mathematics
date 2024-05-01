#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

double lagrangeDerivative(int numPoint, const std::vector<std::pair<double, double>>& tabulatedPoints, double step) {
    int countPoints = tabulatedPoints.size() - 1;
    double result = 0.0;

    for (int i = 0; i <= countPoints; ++i) {
        double pointMult = tabulatedPoints[i].second;

        auto diffMultPart = [&](int a, int b) -> double {
            double subMult = 1.0;
            for (int j = a; j < b; ++j) {
                subMult *= (i - j) * step;
            }
            return subMult;
        };

        double diffMult = diffMultPart(0, i) * diffMultPart(i + 1, countPoints + 1);

        auto gridMultPart = [&](int a, int b) -> double {
            double altMult = 0.0;
            for (int j = a; j < b; ++j) {
                double subMult = 1.0;
                for (int j1 = 0; j1 < std::min(i, j); ++j1) {
                    subMult *= (numPoint - j1);
                }
                for (int j1 = std::min(i, j) + 1; j1 < std::max(i, j); ++j1) {
                    subMult *= (numPoint - j1);
                }
                for (int j1 = std::max(i, j) + 1; j1 <= countPoints; ++j1) {
                    subMult *= (numPoint - j1);
                }
                altMult += subMult;
            }
            return altMult;
        };

        double gridMult = gridMultPart(0, i) + gridMultPart(i + 1, countPoints + 1);

        result += pointMult / diffMult * gridMult;
    }

    return result / step;
}

double mathFunction(double x, int derivativeOrder = 0) {
    if (derivativeOrder == 0) {
        return 0.5 * x * x + cos(2 * x);
    } else if (derivativeOrder == 1) {
        return x - 2 * sin(2 * x);
    }
    return 2 * cos(2 * x + ((derivativeOrder - 2) * M_PI) / 2);
}

std::pair<double, double> residualTerm(int numPoint, double step, int countPoints) {
    double xi = numPoint * step; // Assuming equidistant points
    double factorial = std::tgamma(countPoints + 2); // factorial using gamma function

    double derivativeMax = 0.0;
    for (int i = 0; i <= countPoints; ++i) {
        double derivative = std::fabs(mathFunction(xi, countPoints + 1));
        if (derivative > derivativeMax) {
            derivativeMax = derivative;
        }
    }

    double term = derivativeMax * std::pow(step, countPoints + 1) / factorial;
    return {term, term};
}

int main() {
    int roundingPrecision = 10;
    int derivativeOrder = 1, numberOfPoints = 5, pointIndex = 5;
    std::pair<double, double> interval = {0.6, 1.1};
    double stepSize = (interval.second - interval.first) / numberOfPoints;
    std::vector<std::pair<double, double>> tabulatedPoints;

    for (double point = interval.first; point <= interval.second; point += stepSize) {
        tabulatedPoints.emplace_back(point, mathFunction(point));
    }

    double lagrangeResult = lagrangeDerivative(pointIndex, tabulatedPoints, stepSize);
    double functionResult = mathFunction(tabulatedPoints[pointIndex].first, derivativeOrder);

    auto [minError, maxError] = residualTerm(pointIndex, stepSize, numberOfPoints);

    std::cout << "Lagrange:\t" << std::round(lagrangeResult * std::pow(10, roundingPrecision)) / std::pow(10, roundingPrecision) << std::endl
              << "Function's derivative value:\t" << std::round(functionResult * std::pow(10, roundingPrecision)) / std::pow(10, roundingPrecision) << std::endl
              << "Difference:\t" << std::round(std::abs(lagrangeResult - functionResult) * std::pow(10, roundingPrecision)) / std::pow(10, roundingPrecision) << std::endl
              << std::endl
              << "Minimum error:\t" << minError << std::endl
              << "Maximum error:\t" << maxError << std::endl
              << std::endl
              << "Does the error fall within the range?:\t";
                if((minError < std::abs(lagrangeResult - functionResult) && std::abs(lagrangeResult - functionResult) < maxError)==0)
                    std::cout<<"False";
                else
                    std::cout<<"True";
    return 0;
}
