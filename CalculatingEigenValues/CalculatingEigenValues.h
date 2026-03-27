#pragma once
#include <cmath>
#include <vector>
#include "../External/Matrix/Matrix/Matrix.h"
#include "../External/Matrix/Matrix/Vector.h"
#include "../External/Matrix/Matrix/LinearSystemSolve.h"
using namespace std;

/// <summary>
/// Ищет максимальное собственное число матрицы и соответствующий вектор 
/// с помощью прямого степенного метода
/// </summary>
/// <param name="eigenvalue">Собственное число</param>
/// <param name="x">Собственный вектор</param>
/// <param name="A">Матрица</param>
/// <param name="n">Размерность</param>
/// <param name="tol">Точность</param>
/// <param name="max_iter">Максимальное исло итераций</param>
void powerMethod(double& eigenvalue, Vector& x, Matrix A, int n, double tol, int max_iter = 100000) {
    Vector x_old(n);

    int iter = 0;

    do {
        x_old = x;

        Vector y = A * x;

        double norm = InfVectorNorm(y);
        if (norm < 1e-12) break;

        x = (1 / norm) * y;

        iter++;

    } while (InfVectorNorm(x - x_old) > tol && iter < max_iter);

    Vector Ax = A * x;
    double num = 0.0;
    double den = 0.0;

    for (int i = 0; i < n; i++) {
        num += x.numbersVector[i] * Ax.numbersVector[i];
        den += x.numbersVector[i] * x.numbersVector[i];
    }

    eigenvalue = num / den;
}

/// <summary>
/// Ищет максимальное собственное число матрицы и соответствующий вектор 
/// с помощью обратного степенного метода 
/// </summary>
/// <param name="sigma">Первое приближение</param>
/// <param name="x">Собственный вектор</param>
/// <param name="A">Матрица</param>
/// <param name="n">Размерность</param>
/// <param name="tol">Точность</param>
/// <param name="max_iter">Максимальное исло итераций</param>ы
void inversePowerMethod(double& sigma, Vector& x, Matrix A, int n, double tol, int max_iter = 100000) {
    Matrix I = Matrix::MakeIdentityMatrix(n);

    Vector x_old(n);

    int iter = 0;

    while (iter < max_iter) {
        x_old = x;

        Matrix B = A - sigma * I;

        Vector y = PLU(B, x, n);

        double norm = InfVectorNorm(y);
        if (norm < 1e-12) break;

        x = (1 / norm) * y;

        Vector Ax = A * x;
        double num = 0.0;
        double den = 0.0;

        for (int i = 0; i < n; i++) {
            num += x.numbersVector[i] * Ax.numbersVector[i];
            den += x.numbersVector[i] * x.numbersVector[i];
        }

        double new_sigma = num / den;

        if (abs(new_sigma - sigma) < tol)
            break;

        sigma = new_sigma;
        iter++;
    }
}

/// <summary>
/// Приводит матрицу к форме Хессенберга
/// </summary>
/// <param name="A">Матрица</param>
/// <param name="n">Размерность</param>
/// <returns></returns>
Matrix toHessenberg(Matrix A, int n) {
    Matrix H = A;

    Matrix I = Matrix::MakeIdentityMatrix(n);

    for (int i = 0; i < n - 2; i++) {
        Vector v(n);;

        double s = 0.0;
        for (int j = i + 1; j < n; j++) {
            s += H.numbersMatrix[j][i] * H.numbersMatrix[j][i];
        }

        if (s < 1e-12) continue;

        double sign = (H.numbersMatrix[i + 1][i] >= 0) ? 1.0 : -1.0;
        s = sqrt(s) * sign;

        v.numbersVector[i + 1] = H.numbersMatrix[i + 1][i] - s;

        for (int j = i + 2; j < n; j++) {
            v.numbersVector[j] = H.numbersMatrix[j][i];
        }

        double norm = 0.0;
        for (int j = 0; j < n; j++)
            norm += v.numbersVector[j] * v.numbersVector[j];

        norm = sqrt(norm);
        if (norm < 1e-12) continue;

        for (int j = 0; j < n; j++)
            v.numbersVector[j] /= norm;

        Matrix Hh = I - 2 * (v * v);
        H = Hh * H * Hh;
    }

    return H;
}

/// <summary>
/// Выполняет шаг QR алгоритма
/// </summary>
/// <param name="Q">Матрица Q, ортогональная</param>
/// <param name="R">Матрица R, правая треуголная матрица</param>
/// <param name="n">Размерность</param>
void givensQR(Matrix& Q, Matrix& R, int n) {
    Q = Matrix::MakeIdentityMatrix(n);

    for (int i = 0; i < n - 1; i++) {
        double a = R.numbersMatrix[i][i];
        double b = R.numbersMatrix[i + 1][i];

        if (abs(b) < 1e-12) continue;

        double r = sqrt(a * a + b * b);
        double c = a / r;
        double s = -b / r;

        Matrix G = Matrix::MakeIdentityMatrix(n);

        G.numbersMatrix[i][i] = c;
        G.numbersMatrix[i][i + 1] = -s;
        G.numbersMatrix[i + 1][i] = s;
        G.numbersMatrix[i + 1][i + 1] = c;

        R = G * R;
        Q = Q * G.Transpose();
    }
}

/// <summary>
/// Находит сдвиг Уилкинсона для ускорения сходимость QR алгоритма
/// </summary>
/// <param name="A">Матрица</param>
/// <param name="m">Размерность</param>
/// <returns></returns>
double wilkinsonShift(Matrix& A, int m) {
    double d = (A.numbersMatrix[m - 2][m - 2] - A.numbersMatrix[m - 1][m - 1]) / 2.0;
    double b = A.numbersMatrix[m - 1][m - 2];

    double sign = (d >= 0) ? 1.0 : -1.0;

    return A.numbersMatrix[m - 1][m - 1] -
        sign * b * b / (abs(d) + sqrt(d * d + b * b));
}

/// <summary>
/// Находит все собственные числа матрицы QR алгоритмом
/// </summary>
/// <param name="A">Матрица</param>
/// <param name="n">Размерность</param>
/// <param name="tol">Точность</param>
/// <param name="max_iter">Максимальное число итераций</param>
/// <returns></returns>
vector<double> QR_algorithm(Matrix A, int n, double tol, int max_iter = 100000) {
    Matrix H = toHessenberg(A, n);
    vector<double> eigenvalues;

    int m = n;

    while (m > 1) {
        int iter = 0;

        while (abs(H.numbersMatrix[m - 1][m - 2]) > tol && iter < max_iter) {

            double mu = wilkinsonShift(H, m);

            for (int i = 0; i < m; i++)
                H.numbersMatrix[i][i] -= mu;

            Matrix Q = Matrix::MakeIdentityMatrix(m);
            Matrix R = H.SubMatrix(m);

            givensQR(Q, R, m);

            H = R * Q;

            for (int i = 0; i < m; i++)
                H.numbersMatrix[i][i] += mu;

            iter++;
        }

        eigenvalues.push_back(H.numbersMatrix[m - 1][m - 1]);
        m--;
    }

    eigenvalues.push_back(H.numbersMatrix[0][0]);

    return eigenvalues;
}