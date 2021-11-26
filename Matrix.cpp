#include "Matrix.h"
#include <iostream>
const Matrix&  Matrix::operator=(const Matrix& right){
    _matrix = right._matrix;
    n_size = right.n_size;
    m_size = right.m_size;
    return *this;
}


std::vector<double>& Matrix::operator[](const int index){
    return _matrix[index];
}

const std::vector<double>& Matrix::operator[](const int index) const{
    return _matrix[index];
}


std::ostream& operator<<(std::ostream& os, const Matrix& matrix){
    for (int i = 0; i < matrix.n_size; ++i) {
        for (int j = 0; j < matrix.m_size; ++j)
            os << matrix[i][j] << " ";
        os << std::endl;
    }
    os << std::endl;
    return os;
}



const Matrix operator+(const Matrix& left, const Matrix& right){
    if(left.n_size != right.n_size || left.m_size != right.m_size){
        throw "Wrong summ! Sizes of matrix not equal!";
    }
    Matrix ans(left.n_size, left.m_size);
    for(int i = 0; i < ans.n_size; ++i){
        for(int j = 0; j < ans.m_size; ++j){
            ans[i][j] = left._matrix[i][j] + right._matrix[i][j];
        }
    }
    return ans;
}

void Matrix::SwapRows(size_t i, size_t j) {
    if(i >= get_n() || j >= get_m()) {
        throw std::out_of_range("");
    }
    std::swap(_matrix[i], _matrix[j]);
}

void Matrix::SwapCols(size_t i, size_t j) {
    if(i >=get_n() || j >= get_m()) {
        throw std::out_of_range("");
    }
    for(size_t k = 0; k < get_n(); ++k) {
        std::swap(_matrix[k][i], _matrix[k][j]);
    }
}



const Matrix operator-(const Matrix& left, const Matrix& right){
    if(left.n_size != right.n_size || left.m_size != right.m_size){
        throw "Wrong minus! Sizes of matrix not equal!";
    }
    Matrix ans(left.n_size, left.m_size);
    for(int i = 0; i < ans.n_size; ++i){
        for(int j = 0; j < ans.m_size; ++j){
            ans[i][j] = left._matrix[i][j] - right._matrix[i][j];
        }
    }
    return ans;
}

const Matrix operator*(double left, const Matrix& right){
    Matrix ans = right;
    for(int i = 0; i < ans.n_size; ++i){
        for(int j = 0; j < ans.m_size; ++j){
            ans[i][j] *= left;
        }
    }
    return ans;
}

const Matrix operator*(const Matrix& left, double right){
    return right * left;
}




const Matrix operator*(const Matrix& left, const Matrix& right){
    if(left.m_size != right.n_size){
        throw "Wrong multiply! Sizes of matrix not equal!";
    }
    Matrix ans(left.n_size, right.m_size);
    for(int i = 0; i < ans.n_size; ++i){
        for(int j = 0; j < ans.m_size; ++j){
            for(int k = 0; k < left.m_size; ++k){
                ans[i][j] += left._matrix[i][k] * right._matrix[k][j];
            }
        }
    }
    return ans;
}
const std::vector<double> operator/(const std::vector<double>& v, double val) {
    std::vector<double> res;
    for (int i = 0; i < v.size(); ++i) {
        res[i] = v[i] / val;
    }
    return res;
}

const Matrix operator-(const Matrix& left, double val) {
    Matrix ans(left.n_size, left.m_size);
    for(int i = 0; i < ans.n_size; ++i){
        for(int j = 0; j < ans.m_size; ++j){
            ans[i][j] = left._matrix[i][j] - val;
        }
    }
    return ans;
}


const std::vector<double> operator*(const Matrix& left, const std::vector<double>& right){
    if(left.m_size != (int)right.size()){
        throw "Wrong multiply! Sizes of matrix not equal!";
    }
    std::vector<double> ans(right.size());
    Matrix tmp1(left.get_n(), left.get_m());
    Matrix tmp2(left.get_n(), left.get_m());
    for (int i = 0; i < tmp1.get_n(); ++i) {
        tmp1[i][0] = right[i];
    }
    tmp2 = left * tmp1;
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = tmp2[i][0];
    }
    return ans;
}



int Matrix::get_m() const{
    return m_size;
}

int Matrix::get_n() const{
    return n_size;
}

void Matrix::make_ones(){
    if(!is_quadratic()){
        throw "Ones can be done only with quadratic matrix";
    }
    _matrix.assign(n_size, std::vector<double>(m_size, 0.0));
    for(int i = 0; i < n_size; ++i){
        _matrix[i][i] = 1.0;
    }
}

void Matrix::transpose(){
    std::vector<std::vector<double>> temp(m_size, std::vector<double>(n_size));
    for(int i = 0; i < n_size; ++i){
        for(int j = 0; j < m_size; ++j){
            temp[j][i] = _matrix[i][j];
        }
    }
    std::swap(n_size, m_size);
    _matrix.swap(temp);
}

// Create 1 x 1 matrix with 0.0 element
Matrix::Matrix(){
    _matrix.assign(1, std::vector<double>(1, 0));
    n_size = m_size = 1;
}

// Create n x m matrix with 0.0 elements
Matrix::Matrix(int n, int m){
    if(!n || !m){
        throw "Matrix size must be > 0";
    }
    _matrix.assign(n, std::vector<double>(m, 0));
    n_size = n;
    m_size = m;
}


bool Matrix::is_quadratic() const{
    return n_size == m_size;
}


bool Matrix::is_three_diagonal() const{
    if(!is_quadratic()){
        return false;
    }
    for(int i = 0; i < n_size; ++i){
        for(int j = 0; j < m_size; ++j){
            if((abs(i - j) > 1) && _matrix[i][j]){
                return false;
            }
        }
    }
    return true;
}

bool Matrix::is_simmetric() const{
    if(!is_quadratic()){
        return false;
    }
    for(int i = 0; i < n_size; ++i){
        for(int j = i + 1; j < m_size; ++j){
            if(_matrix[i][j] != _matrix[j][i]){
                return false;
            }
        }
    }
    return true;
}

bool FindParityOfPermutation(Matrix &P) {
    int permutation_count = 0;
    for(int i = 0; i < P.get_n(); ++i) {
        if(P.GetElement(i,i) != 1) {
            ++permutation_count;
        }
    }
    return permutation_count % 2;
}

void Matrix::SetElement(int i, int j, double value) {
    if(i >= get_n() || j >= get_m()) {
        throw std::out_of_range("");
    }
    _matrix[i][j] = value;
}

double &Matrix::GetElement(int i, int j) {
    if(i >= get_n() || j >= get_m()) {
        throw std::out_of_range("");
    }
    return _matrix[i][j];
}

const double &Matrix::GetElement(int i, int j) const {
    if(i >= get_n() || j >= get_m()) {
        throw std::out_of_range("");
    }
    return _matrix[i][j];
}



double Matrix::get_norm() const{
    double max = 0.0;
    for(int i = 0; i < n_size; ++i){
        double ans = 0.0;
        for(int j = 0; j < m_size; ++j){
            ans += abs(_matrix[i][j]);
        }
        max = max > ans ? max : ans;
    }
    return max;
}

void Matrix::SetCol(const Matrix &vec, int col_number) {
    if (vec.get_m() == 1) {
        for (int i = 0; i < get_n(); ++i) {
            _matrix[i][col_number] = vec[i][0];
        }
    }
}
Matrix Matrix::GetCol(int col_number) {
    Matrix result(get_n(), 1);
    for (int i = 0; i < get_n(); ++i) {
        result[i][0] = GetElement(i, col_number);
    }
    return result;
}
double Matrix::get_upper_norm() const{
    double max = 0.0;
    for(int i = 0; i < n_size; ++i){
        double ans = 0.0;
        for(int j = 0; j <= i; ++j){
            ans += abs(_matrix[i][j]);
        }
        max = max > ans ? max : ans;
    }
    return max;
}

double scal_prod(const Matrix& v1, const Matrix& v2) {
    double ans;
    if (v1.get_m() == 1 && v2.get_m() == 1) {
        for (int i = 0; i < v1.get_n(); ++i) {
            ans += v1[i][0] * v2[i][0];
        }
        return ans;
    }
}

double VectNorm(const Matrix& v) {
    if(v.get_m() == 1 ) {
        double res;
        for(int i = 0; i < v.get_n(); ++i) {
            res += v[i][0] * v[i][0];
        }
        return sqrt(res);
    }
}


