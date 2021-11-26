#include "Matrix.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#define value_pair pair<pair<double, double>, pair<double, double>>
using namespace std;
void print_vector_x(const vector<pair<double, double>>& x){
    for(unsigned i = 0; i < x.size(); ++i){
        cout << 'l' << i + 1 << " = " << x[i].first;
        if(x[i].second){
            if(x[i].second > 0){
                cout << " + ";
            }else{
                cout << " - ";
            }
            cout << abs(x[i].second) << "i";
        }
        cout << endl;
    }
    cout << endl;
}

void print_solution(const vector<pair<double, double>>& x, int itter){
    cout << "Values:" << endl;
    print_vector_x(x);
    cout << "Itterations: " << itter << endl;
}

double mult_1xn_nx1_vecs(const vector<double>& left, const vector<double>& right){
    double ans = 0.0;
    if(left.size() != right.size()){
        throw "Wrong sizes of vectors!";
    }
    for(unsigned i = 0; i < left.size(); ++i){
        ans += right[i] * left[i];
    }
    return ans;
}

Matrix mult_nx1_1xn_vecs(const vector<double>& left, const vector<double>& right){
    Matrix ans(left.size(), right.size());
    for(int i = 0; i < ans.get_n(); ++i){
        for(int j = 0; j < ans.get_m(); ++j){
            ans[i][j] += left[i] * right[j];
        }
    }
    return ans;
}

double sign(double num){
    if(!num){
        return 0.0;
    }
    return num > 0 ? 1.0 : -1.0;
}

void solve_eq(double a, double b, double c, pair<pair<double, double>, pair<double, double>>& ans){
    double D = b * b - 4.0 * a * c;
    if(D >= 0.0){
        ans.first.first = (-b + sqrt(D)) / (2 * a);
        ans.first.second = 0.0;
        ans.second.first = (-b - sqrt(D)) / (2 * a);
        ans.second.second = 0.0;
    }else{
        ans.first.first = -b / (2 * a);
        ans.first.second = sqrt(-D) / (2 * a);
        ans.second.first = -b / (2 * a);
        ans.second.second = -sqrt(-D) / (2 * a);;
    }
}

double complex_check(const value_pair& last, const value_pair& cur){
    pair<double, double> r1, r2;
    // x[j]
    r1.first = cur.first.first - last.first.first; // x
    r1.second = cur.first.second - last.first.second; // y
    // x[j+1]
    r2.first = cur.second.first - last.second.first; // x
    r2.second = cur.second.second - last.second.second; // y
    return max(sqrt(r1.first*r1.first + r1.second*r1.second), sqrt(r2.first*r2.first + r2.second*r2.second));
}

void QR_dec(const Matrix& A, Matrix& Q, Matrix& R){
    Matrix E(A.get_n(), A.get_m());
    E.make_ones();
    Q = E;
    R = A;

    for(int j = 0; j < R.get_m() - 1; ++j){
        vector<double> v(R.get_n(), 0.0);
        double norm = 0.0;

        v[j] = R[j][j];
        for(int i = j; i < R.get_n(); ++i){
            norm += R[i][j] * R[i][j];
        }
        norm = sqrt(norm);
        v[j] += sign(R[j][j]) * norm;

        for(int i = j + 1; i < R.get_n(); ++i){
            v[i] = R[i][j];
        }

        // Hausdorf matrix:
        Matrix H = E - (2.0 / mult_1xn_nx1_vecs(v, v)) * mult_nx1_1xn_vecs(v, v);

        Q = Q * H;
        R = H * R;
    }
}

int QRmethod_values(const Matrix& A, vector<pair<double, double>>& x, double alfa){
    Matrix Q, R, A_k = A;
    x.resize(A_k.get_m());
    int itter = 0;
    double check;
    value_pair curr;
    bool flag = true;

    for(itter = 0; flag; ++itter){
        QR_dec(A_k, Q, R); // decompose
        A_k = R * Q;

        flag = false;
        for(int j = 0; j < A_k.get_m(); ++j){
            check = 0.0;
            for(int i = j + 1; i < A_k.get_n(); ++i){
                check += A_k[i][j] * A_k[i][j];
            }
            check = sqrt(check);

            if(check > alfa){
                solve_eq(1.0, -A_k[j][j]-A_k[j+1][j+1], A_k[j][j]*A_k[j+1][j+1]-A_k[j+1][j]*A_k[j][j+1], curr);
                if(complex_check(curr, value_pair(x[j], x[j+1])) > alfa){
                    flag = true;
                }
                x[j] = curr.first;
                x[j + 1] = curr.second;
                ++j;
            }else{
                x[j].first = A_k[j][j];
                x[j].second = 0.0;
            }
        }
    }
    return itter;
}

int main() {
    fstream is("QR.txt");
    if (is) {
        Matrix A, Q, R;
        vector<pair<double, double>> x;
        double accuracy = 0.01;
        int size;
        is >> size;
        cout << size << endl;
        A = Matrix(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                is >> A[i][j];
            }
        }
        cout << "Matrix A:" << endl;
        cout << A << endl;
        QR_dec(A, Q, R);
        cout << Q << endl;
        cout << R << endl;
        cout << "Q*R:" << endl;
        cout << Q * R << endl;
        cout << "QR Method:" << endl;
        int itter = QRmethod_values(A, x, accuracy);
        print_solution(x, itter);
    } else {
        cout << "Bad source!" << endl;
    }
    return 0;
}