#include <vector>

using namespace std;

void givens_rotation(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N);

vector <double> reverse_Gaussian(vector <vector <double>> & matrix_a, vector <double> & phi, int N);

void scan_m(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N);

void print_m(vector <double> & x, int N);
