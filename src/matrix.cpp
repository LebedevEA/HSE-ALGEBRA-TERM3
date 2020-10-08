#include <memory>
#include <cmath>
#include "matrix.h"

Matrix::Matrix()
    : rows_(0)
    , columns_(0)
    , table_(nullptr)
{}

Matrix::Matrix(const Matrix& other)
    : rows_(other.rows_)
    , columns_(other.columns_)
    , table_(nullptr) {
    if (!(other.rows_ == 0 or other.columns_ == 0)) {
        table_ = std::make_unique<double*[]>(other.rows_);
        std::unique_ptr<double[]> tmp = std::make_unique<double[]>(other.rows_ * other.columns_);
        for (int i = 0; i < this->rows_; i++)
            table_[i] = &tmp[this->columns_ * i];
        for (int c = 0; c < this->columns_; c++)
            for (int r = 0; r < this->rows_; r++)
                table_[r][c] = other[r][c];
        tmp.release();
    }
}

Matrix::Matrix(int rows, int columns)
    : rows_(rows)
    , columns_(columns) {
    if (!(rows == 0 or columns_ == 0)) {
        table_ = std::make_unique<double*[]>(rows_);
        std::unique_ptr<double[]> tmp = std::make_unique<double[]>(rows_ * columns_);
        for (int i = 0; i < this->rows_; i++)
            table_[i] = &tmp[this->columns_ * i];
        for (int c = 0; c < this->columns_; c++)
            for (int r = 0; r < this->rows_; r++)
                table_[r][c] = 0;
        tmp.release();
    }
}


Matrix::~Matrix() {
    if (table_) delete[] table_[0];
}

Matrix& Matrix::operator=(Matrix other) {
    this->swap(other);
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& other) {
    Matrix m = *this + other;
    *this = m;
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& other) {
    Matrix m = *this * other;
    *this = m;
    return *this;
}

Matrix operator+(const Matrix& fst, const Matrix& snd) {
    if (fst.rows_ != snd.rows_ or fst.columns_ != snd.columns_)
        throw MatrixException("ADD: dimensions do not match.\n");
    Matrix m(fst);
    for (int c = 0; c < fst.columns_; c++)
        for (int r = 0; r < fst.rows_; r++)
            m[r][c] += snd[r][c];
    return m;
}

Matrix::MatrixRow Matrix::operator[](int row) const {
    if (row >= rows_ or row < 0 or columns_ == 0 or rows_ == 0)
        throw MatrixException("ACCESS: bad index.\n");
    return MatrixRow(columns_, table_[row]);
}

Matrix operator*(const Matrix& fst, const Matrix& snd) {
    if (fst.columns_ != snd.rows_)
        throw MatrixException("MUL: #arg1.columns != #arg2.rows.\n");
    Matrix m(fst.rows_, snd.columns_ );
    for (int i = 0; i < fst.rows_; i++)
        for (int j = 0; j < snd.columns_; j++)
            for (int k = 0; k < snd.rows_; k++)
                m[i][j] += fst[i][k] * snd[k][j];
    return m;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (int i = 0; i < m.rows_; i++) {
        for (int j = 0; j < m.columns_ - 1; j++)
            os  << std::setw(7) << m[i][j] << ' ';
        os << std::setw(7) << m[i][m.columns_ - 1] << '\n';
    }
    return os;
}

std::ifstream& operator>>(std::ifstream& is, Matrix& m) {
    if (!is.is_open())
        throw MatrixException("LOAD: unable to open file.\n");
    int rows = 0, columns = 0;
    int counter = 0;
    is >> rows >> columns;
    if (rows == 0 or columns == 0)
        throw MatrixException("LOAD: invalid file format.\n");
    Matrix mtx(rows, columns);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++) {
            if (is >> mtx[i][j]) counter++;
        }
    if (counter != rows * columns)
        throw MatrixException("LOAD: invalid file format.\n");
    m = mtx;
    return is;
}

void Matrix::swap(Matrix& other) noexcept {
    std::swap(this->rows_, other.rows_);
    std::swap(this->columns_, other.columns_);
    std::swap(this->table_, other.table_);
}

Matrix::MatrixColumn Matrix::getColumn(int num) {
    return MatrixColumn(*this, num);
}

int Matrix::rows() const {
    return rows_;
}

int Matrix::columns() const {
    return columns_;
}

Matrix::Matrix(MatrixType type, int size)
    : Matrix(size, size) {
    if (type == ID) {
        for (int i = 0; i < size; i++) {
            (*this)[i][i] = 1.0;
        }
    }
}

Matrix::Matrix(Matrix&& mtx)
    : Matrix() {
    swap(mtx);
}

Matrix Matrix::transpose() const {
    auto A = Matrix(this->columns(), this->rows());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.columns(); j++) {
            A[i][j] = (*this)[j][i];
        }
    }
    return A;
}

Matrix::Matrix(const Matrix::MatrixColumn& mtx) : Matrix(mtx.rows(), 1) {
    for (int row = 0; row < rows(); row++) {
        (*this)[row][0] = mtx[row];
    }
}

Matrix::MatrixRow::MatrixRow(int max_rows, double* column) noexcept
    : max_columns_(max_rows)
    , row_(column)
{}

double& Matrix::MatrixRow::operator[](int column) {
    if (column >= max_columns_ or column < 0)
        throw MatrixException("ACCESS: bad index.\n");
    return row_[column];
}

Matrix& MatrixArray::operator[](std::size_t index) {
    if (index >= size_)
        throw MatrixException("ACCESS: bad index.\n");
    return ma_[index];
}

MatrixArray::MatrixArray(std::size_t size) : size_(size) , ma_(new Matrix[size]) {}

MatrixArray::~MatrixArray() {
    delete[] ma_;
}

Matrix::MatrixColumn::MatrixColumn(Matrix& ma, int col_num) noexcept
    : ma_(ma)
    , column_(col_num)
{}

double& Matrix::MatrixColumn::operator[](int row) const {
    return ma_[row][column_];
}

int Matrix::MatrixColumn::rows() const {
    return ma_.rows();
}

double product(const Matrix::MatrixColumn& first, const Matrix::MatrixColumn& second) {
    double ret_val = 0;
    if (first.rows() != second.rows()) {
        throw MatrixException("Product: dimensions does not match.");
    }
    for (int i = 0; i < first.rows(); i++) {
        ret_val += first[i] * second[i];
    }
    return ret_val;
}


double norm(const Matrix::MatrixColumn& col) {
    return std::sqrt(product(col, col));
}