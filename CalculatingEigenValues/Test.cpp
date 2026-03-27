#include <iostream>
#include <cmath>
#include <locale>
#include "CalculatingEigenValues.h"

using namespace std;

int main() {
    double tol = 1e-8;
    setlocale(LC_ALL, "Russian");

    while (true) {
        int n;
        cout << "Введите размерность матрицы (0 - выход): ";
        cin >> n;

        if (n == 0) break;

        Matrix A(n);
        cout << "Введите матрицу:\n";
        A.InputFromConsole();

        while (true) {
            cout << "\nВыберите метод:\n";
            cout << "1 - Степенной метод\n";
            cout << "2 - Обратный степенной метод\n";
            cout << "3 - QR алгоритм (все собственные значения)\n";
            cout << "4 - Ввести новую матрицу\n";
            cout << "0 - Выход\n";

            int choice;
            cin >> choice;

            if (choice == 0) return 0;   

            if (choice == 4) break;

            if (choice == 1) {
                Vector x(n);
                vector<double> init(n, 1.0);
                x.numbersVector = init;

                double eigen = 0.0;

                powerMethod(eigen, x, A, n, tol);

                cout << "Максимальное eigenvalue: " << eigen << endl;
                cout << "Собственный вектор:\n";
                double norma = L2VectorNorm(x);
                x = (1 / norma) * x;
                x.Print();
            }

            else if (choice == 2) {
                double sigma;
                cout << "Введите начальное приближение sigma: ";
                cin >> sigma;

                Vector x(n);
                vector<double> init(n, 1.0);
                x.numbersVector = init;

                inversePowerMethod(sigma, x, A, n, tol);

                cout << "Найденное eigenvalue: " << sigma << endl;
                cout << "Собственный вектор:\n";
                double norma = L2VectorNorm(x);
                x = (1 / norma) * x;
                x.Print();
            }

            else if (choice == 3) {
                vector<double> eigenvalues = QR_algorithm(A, n, tol);

                cout << "Собственные значения:\n";
                for (double val : eigenvalues)
                    cout << val << " ";
                cout << endl;
            }

            else {
                cout << "Неверный выбор\n";
            }
        }
    }

    return 0;
}