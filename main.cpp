#include <iostream>
#include <fstream>
#include "givens_rotation.h"

using namespace std;

int main()
{
    vector <double> x, phi_qr, phi_ne; 
    vector <vector <double>> matrix_a_qr, matrix_a_ne;

    fstream file;
    file.open("experiment_3.txt");
    while (!file.eof())
    {
        double x_i, phi_i;
        file >> x_i >> phi_i;
        x.push_back(x_i);
        phi_qr.push_back(phi_i);
        phi_ne.push_back(phi_i);
    }
    file.close(); 

    int N = 5;

    matrix_a_qr = fill_matrix_a(x, N);
    matrix_a_ne = fill_matrix_a(x, N);

    print_m(matrix_a_qr);
    print_m(matrix_a_ne);

    givens_rotation(matrix_a_qr, phi_qr, N); // переопределили A и phi
    
    vector <double> c_qr, c_ne;

    c_qr = reverse_Gaussian(matrix_a_qr, phi_qr, N);
    c_ne = normal_equation(matrix_a_ne, phi_ne, N);
    
    double h = x.back() / 1000.; // шаг сетки
    vector <double> y_qr(1000), y_ne(1000), t(1000);
    t[0] = 0;
    y_qr[0] = f(t[0], c_qr, N);
    y_ne[0] = f(t[0], c_ne, N);
    for (int i = 1; i < 1000; ++i)
    {
        t[i] = t[i - 1] + h;
        y_qr[i] = f(t[i], c_qr, N);
        y_ne[i] = f(t[i], c_ne, N);
    }
    
    ofstream file_r("result.txt");
    for (int i = 0; i < 1000; ++i)
        file_r <<  t[i] << ' ' << y_ne[i] << ' ' << y_qr[i] << endl;

    file_r.close();
    

    return 0;
}