#include <iostream>
#include <vector>
#include <cmath>
#include "givens_rotation.h"

using namespace std;

void givens_rotation(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N) 
{
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
    double buffer_1, buffer_2;
    vector <double> matrix_x(N);
    
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

void scan_m(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N)
{
    cout << "Enter the matrix A: " << endl;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < N; ++j)
        {
            cin >> matrix_a[i][j];
        }
  
    cout << "Enter the matrix b: " << endl;
    for (int i = 0; i < m; ++i)
    {
        cin >> phi[i];
    }   
}
void print_m(vector <double> & x, int N) 
{
   cout << "Result: " << endl;
    for (int i = 0; i < N; ++i)
        cout << x[i] << endl;
}

