#ifndef NM_MATRIX_H
#define NM_MATRIX_H


#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


class Matrix{
public:
    Matrix();
    Matrix(int n, int m);

    void SwapRows(size_t i, size_t j);
    void SwapCols(size_t i, size_t j);
    void make_ones();
    void transpose();
    void SetElement(int i, int j, double value);
    void GetElement(int i, int j, double value);
    double &GetElement(int i, int j);
    const double &GetElement(int i, int j) const;
    void SetCol(const Matrix& vec, int col_number);
    Matrix GetCol(int col_number);
    std::vector<double>& operator[](const int index);
    const std::vector<double>& operator[](const int index) const;

    friend const Matrix operator+(const Matrix& left, const Matrix& right);
    friend const Matrix operator-(const Matrix& left, const Matrix& right);
    friend const Matrix operator-(const Matrix& left, double val);
    friend const Matrix operator*(const Matrix& left, const Matrix& right);
    friend const std::vector<double> operator*(const double val, const std::vector<double>& right);
    friend const std::vector<double> operator*(const Matrix& left, const std::vector<double>& right);
    friend const std::vector<double> operator/(const std::vector<double>& v, double val);
    friend const Matrix operator*(double left, const Matrix& right);

    const Matrix& operator=(const Matrix& right);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);



    double get_norm() const;
    double get_upper_norm() const;

    int get_n() const;
    int get_m() const;

    bool is_three_diagonal() const;
    bool is_simmetric() const;

    bool is_quadratic() const;



private:
    std::vector<std::vector<double>> _matrix;
    int n_size;
    int m_size;
};
bool FindParityOfPermutation(Matrix &P);
double scal_prod(const Matrix& v1, const Matrix& v2);
double VectNorm(const Matrix& m);
#endif
