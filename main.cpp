#include <iostream>
#include "Matrix.h"
void arnoldi(const Matrix& A, const Matrix& v,
             size_t k, Matrix& V, Matrix& H) {
    int n = A.get_m();
    Matrix vt(n, 1);
    V.SetCol((1 / VectNorm(v)) * v , 0  );
    for (int m = 0; m < k; m++) {
        vt = A * V.GetCol(m);
        for (int j=0; j < m+1; j++) {
            H[j][m] = scal_prod(vt, (V.GetCol(j)));
            vt = vt - H[j][m] * V.GetCol(j) ;
        }
        H[m+1][m] = VectNorm(vt);
        double eps = 1e-5;
        if (H[m+1][m] > eps) {
            V.SetCol(1 / H[m+1][m] *  vt, m+1);
        } else {
            break;
        }
    }
}
int main() {
    std::fstream is("matrix.txt");
    if (is) {
        Matrix A;
        int k = 9;
        int size;
        is >> size;
        std::cout << size << std::endl;
        A = Matrix(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                is >> A[i][j];
            }
        }
        std::cout << A << std::endl;
        Matrix vect(size, 1);
        for (int i = 0; i < size; ++i) {
            is >> vect[i][0];
        }
        std::cout << vect << std::endl;
        Matrix V(size, k);
        Matrix H(k + 1, k );
        arnoldi(A, vect, k, V, H);
        std::cout << V << std::endl;
        std::cout << H << std::endl;
    }
}