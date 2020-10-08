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

std::vector<std::pair<Matrix, double>> algoQR(Matrix A, int itNo = 5000) {
    Matrix limQ(ID, A.rows()), INIT = A;
    while (itNo --> 0) {
        auto [Q, R] = decomposeQR(A);
        A = R * Q;
        limQ *= Q;
    }
    std::cout << limQ;
    std::vector<std::pair<Matrix, double>> ret_vec;
    for (int i = 0; i < limQ.columns(); i++) {
        Matrix v = Matrix(limQ.getColumn(i));
        auto u = INIT * v;
        int cnt = 0;
        double lambda = 0;
        // вообще поиск чисел можно было заменить
        // на взятие чисел с диагонали матрицы А
        // но это быстрее пришло в голову
        for (int j = 0; j < u.rows(); j++) {
            if (u[j][0] > 1e-10) {
                lambda += u[j][0] / v[j][0];
                cnt++;
            }
        }
        lambda /= cnt;
        ret_vec.emplace_back(v, lambda);
    }
    return ret_vec;
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
         * Если название файла не mtx.txt, то раскоментировать 109 строку
         * и потом после запуска вписать название
         */
        std::cout << std::fixed << std::setprecision(3) << std::setfill(' ');
        std::cout << "Enter file name: ";
        std::string filename("mtx.txt");
//        std::cin >> filename;
        std::cout << std::endl;
        std::ifstream fin(filename);
        Matrix A;
        fin >> A;
        auto vec = algoQR(A);
        for (auto [v, lambda] : vec) {
            std::cout << "lambda = " << lambda << '\n' << v << /*'\n' << A * v << */"\n\n";
        }
    } catch (std::exception& e) {
        std::cout << e.what();
    }
    return 0;
}
