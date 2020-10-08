#pragma once

#include <iostream>
#include <fstream>
#include <memory>

using std::swap;

enum MatrixType {
    ID
};

class Matrix {
public:
    class MatrixColumn {
    public:
        MatrixColumn(Matrix& ma, int col_num) noexcept;
        double& operator[](int row) const;
        [[nodiscard]] int rows() const;
    private:
        const Matrix& ma_;
        int column_;
    };
    class MatrixRow {
    public:
        MatrixRow(int max_rows, double* column) noexcept;
        double& operator[](int row);
    private:
        int max_columns_;
        double* row_;
    };
    Matrix();
    Matrix(const MatrixColumn& mtx);
    Matrix(const Matrix&);
    Matrix(Matrix&& mtx);
    Matrix(int columns, int rows);
    Matrix(MatrixType type, int size);
    ~Matrix();
    [[nodiscard]] Matrix transpose() const;
    [[nodiscard]] MatrixColumn getColumn(int num);
    [[nodiscard]] int rows() const;
    [[nodiscard]] int columns() const;
    Matrix& operator=(Matrix other);
    Matrix& operator+=(const Matrix& other);
    Matrix& operator*=(const Matrix& other);
    MatrixRow operator[](int column) const;
private:
    int rows_, columns_;
    std::unique_ptr<double*[]> table_;

    void swap(Matrix& other) noexcept;

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend std::ifstream& operator>>(std::ifstream& is, Matrix& m);
    friend Matrix operator+(const Matrix& fst, const Matrix& snd);
    friend Matrix operator*(const Matrix& fst, const Matrix& snd);
};

class MatrixException : public std::logic_error { // Другого с конструтором от char* не нашел (ну, нашел, но названия там были не лучше)
public:
    using std::logic_error::logic_error;
};

class MatrixArray {
public:
    MatrixArray() = delete;
    explicit MatrixArray(std::size_t size);
    ~MatrixArray();
    MatrixArray& operator=(const MatrixArray&) = delete;
    Matrix& operator[](std::size_t index);
private:
    std::size_t size_;
    Matrix* ma_;
};

double product(const Matrix::MatrixColumn& first, const Matrix::MatrixColumn& second);

double norm(const Matrix::MatrixColumn& col);