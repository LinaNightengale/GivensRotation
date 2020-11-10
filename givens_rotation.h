#include <vector>

using namespace std;

vector <double> polynomyal(double x, int N);

vector <vector <double>> fill_matrix_a(vector <double> & x, int N);

void givens_rotation(vector <vector <double>> & matrix_a, vector <double> & phi, int N);

vector <double> reverse_Gaussian(vector <vector <double>> & matrix_a, vector <double> & phi, int N);

vector <vector <double>> transpose(vector <vector <double>> & matrix_a, int N);

vector <vector <double>> mult(vector <vector <double>> matrix_left, vector< vector<double>> matrix_right);

vector <double> normal_equation(vector <vector <double>> & matrix_a, vector <double> & phi, int N);

double f(double x, vector <double> & c, int N);

void print_m(vector <vector <double>> & m);
