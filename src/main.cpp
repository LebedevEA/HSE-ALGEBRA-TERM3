#include <iostream>
#include <vector>
#include <iomanip>
#include "matrix.h"

struct QR {
    Matrix Q, R;
};

// Я там где-то что-то с индексами намудрил, там ничего интересного нет, можно не разбираться

Matrix mul1(const std::vector<Matrix>& vec) {
    int size = vec.back().columns();
    Matrix M(ID, size);
    for (const Matrix& mtx : vec) {
        for (int row = 0; row < size; row++)
            M[row][size - 1] = mtx[row][size - 1];
        size--;
    }
    return M;
}

Matrix mul3(const std::vector<Matrix>& vec) {
    int size = vec.back().columns();
    Matrix M(ID, size);
    size--;
    while (size >= 0) {
        M[size][size] = vec[vec.size() - size - 1][size][size];
        size--;
    }
    return M;
}

QR decomposeQR(Matrix A) {
    std::vector<Matrix> elementary1, elementary3;
    for (int i = 1; i < A.columns(); i++) {
        Matrix S(ID, A.columns());
        Matrix::MatrixColumn cur(A, i);
        for (int j = 0; j < i; j++) {
            double prod = product(cur, A.getColumn(j)) / product(A.getColumn(j), A.getColumn(j));
            S[j][i] = -prod;
            for (int k = 0; k < A.rows(); k++)
                A[k][i] += S[j][i] * A[k][j];
        }
        for (int j = 0; j < i; j++)
            S[j][i] = - S[j][i];
        elementary1.push_back(std::move(S));
    }
    for (int i = 0; i < A.columns(); i++) {
        Matrix S(ID, A.columns());
        Matrix::MatrixColumn cur(A, i);
        double curNorm = norm(cur);
        S[i][i] = 1 / curNorm;
        for (int j = 0; j < A.rows(); j++)
            A[j][i] *= S[i][i];
        S[i][i] = curNorm;
        elementary3.push_back(std::move(S));
    }
    std::reverse(elementary3.begin(), elementary3.end());
    std::reverse(elementary1.begin(), elementary1.end());
    return { A, mul3(elementary3) * mul1(elementary1) };
}

int main() {
    try {
        /*
         * Формат файлика, пример:
         * 3 2
         * 1 0
         * 2 1
         * 1 4
         * (звездочки тут для красоты, а не для файлика)
         * Первые 2 строки - размерность, дальше сама матрица
         */
        std::cout << std::fixed << std::setprecision(3) << std::setfill(' ');
        std::cout << "Enter file name: ";
        std::string filename("mtx.txt");
//        std::cin >> filename;
        std::cout << std::endl;
        std::ifstream fin(filename);
        Matrix A;
        fin >> A;
        auto [Q, R] = decomposeQR(A);
        std::cout << "A:\n" << A << "\n"
                  << "Q:\n" << Q << "\n"
                  << "R:\n" << R << "\n"
                  << "Q * R:\n" << Q * R << "\n"
                  << "Q * (Q ^ T):\n" << Q * Q.transpose() << "\n";
    } catch (std::exception& e) {
        std::cout << e.what();
    }
    return 0;
}
