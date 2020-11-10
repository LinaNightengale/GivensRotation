#include <iostream>
#include <vector>
#include <cmath>
#include "givens_rotation.h"

using namespace std;

vector <double> polynomyal(double x, int N) // полиномы Чебышева 1го рода
{
    vector <double> p_j = { 1, x };

    for (int j = 2; j < N; ++j)
        p_j.push_back(2. * x * p_j[j - 1] - p_j[j - 2]);

    return p_j;
}

vector <vector <double>> fill_matrix_a(vector <double> & x, int N) // заполняем матрицу А
{
    vector <vector <double>> matrix;

    int m = x.size();

    for (int i = 0; i < m; ++i)
        matrix.push_back(polynomyal(x[i], N));

    return matrix;
}

void givens_rotation(vector <vector <double>> & matrix_a, vector <double> & phi, int N) // вращение Гивенса
{
    int m = matrix_a.size();

    for (int i = 0; i < N; ++i)
    {
        for (int j = (m - 1); j >= 0; --j)
        {
            if (j > i)
            {
                double c = matrix_a[j - 1][i] / sqrt(powf(matrix_a[j - 1][i], 2.) + powf(matrix_a[j][i], 2.));
                double s = -matrix_a[j][i] / sqrt(powf(matrix_a[j - 1][i], 2.) + powf(matrix_a[j][i], 2.));

                vector <double> buffer_1, buffer_2;

                for (int k = 0; k < N; ++k)
                {
                    buffer_1.push_back(c * matrix_a[j - 1][k] - s * matrix_a[j][k]);
                    buffer_2.push_back(s * matrix_a[j - 1][k] + c * matrix_a[j][k]);
                }

                double buffer_3 = c * phi[j - 1] - s * phi[j];
                double buffer_4 = s * phi[j - 1] + c * phi[j];

                matrix_a[j - 1] = buffer_1;
                matrix_a[j] = buffer_2;

                phi[j - 1] = buffer_3;
                phi[j] = buffer_4;
            }
            else
                continue;
        }
    }
}

vector <double> reverse_Gaussian(vector <vector <double>> & matrix_a, vector <double> & phi, int N)
{
    vector <double> matrix_x(N);

    double buffer_1, buffer_2;

    for (int i = N - 1; i >= 0; --i)
    {
        buffer_1 = 0.;
        for (int j = i + 1; j < N; ++j)
        {
            buffer_2 = matrix_a[i][j] * matrix_x[j];
            buffer_1 += buffer_2;
        }
        matrix_x[i] = (phi[i] - buffer_1) / matrix_a[i][i];
    }
    return matrix_x;
}

vector <vector <double>> transpose(vector <vector <double>> & matrix_a, int N)
{
    int m = matrix_a.size();

    vector <vector <double>> matrix_a_t(N, vector <double> (m, 0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j)
            matrix_a_t[i][j] = matrix_a[j][i];

    return matrix_a_t;
}

vector <vector <double>> mult(vector <vector <double>> matrix_left, vector< vector<double>> matrix_right)
{
  int l = matrix_left.size(), m = matrix_left[0].size(), n = matrix_right[0].size();

  vector <vector <double>> matrix(l, vector <double> (n, 0));

  for (int i = 0; i < l; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < m; ++k)
        matrix[i][j] += matrix_left[i][k] * matrix_right[k][j];

  return matrix;
}

vector <double> normal_equation(vector <vector <double>> & matrix_a, vector <double> & phi, int N)
{
    int m = phi.size();
    vector <double> x;
    vector <vector <double>> AT = transpose(matrix_a, N);
    vector <vector <double>> ATA = mult(AT, matrix_a);
    vector <double> ATphi (N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j)
        {
            ATphi[i] += AT[i][j] * phi[j];
        }
    givens_rotation(ATA, ATphi, N);
    x = reverse_Gaussian(ATA, ATphi, N);
    return x;
}

double f(double x, vector <double> & c, int N)
{
    double y = 0.;
    vector <double> buffer = polynomyal(x, N);

    for (int j = 0; j < N; ++j)
    {
        y += c[j] * buffer[j];
    }
    return y;
}

void print_m(vector <vector <double>> & m) 
{
    for (auto col : m) 
    {
        for (auto el : col)
            cout << el << ' ';
        cout << endl;
    }
}