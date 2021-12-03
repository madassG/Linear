#include <math.h>

#include "linear.h"
#include "exceptions.cpp"
#include "operators.cpp"

template <typename T>
Linear::Matrix<T>::Matrix(const std::vector<std::vector<T>>& matrix) : _matrix(matrix.begin(), matrix.end()),
                                                                       _rows(matrix.size()),
                                                                       _columns(0){
    if (_matrix.empty()) return;

    _columns = _matrix[0].size();
    for (const auto& row : _matrix) {
        if (row.size() != _columns) throw Linear::IncorrectInitialVector("Incorrect initial vector, "
                                                                "row sizes differ from each other");
    }
}

template <typename T>
Linear::Matrix<T>::Matrix(const std::vector<T>& row) : _matrix(1, std::vector<T>(row.begin(), row.end())), _rows(1), _columns(row.size()) {}

template <typename T>
Linear::Matrix<T>::Matrix(int rows, int cols) : _matrix(cols, std::vector(rows, T(0))), _rows(rows), _columns(cols) {}

template <typename T>
Linear::Matrix<T>::Matrix(int N) : _matrix(N, std::vector(N, T(0))), _rows(N), _columns(N) {}

template <typename T>
Linear::Matrix<T>::Matrix() : _matrix(1, std::vector(1, T(0))), _rows(1), _columns(1) {}


template <typename T>
T& Linear::Matrix<T>::At(size_t i, size_t j) {
    if (i >= _rows) throw OutOfRange("Row value " + std::to_string(i) + " is out of range.");
    if (j >= _columns) throw OutOfRange("Column value " + std::to_string(j) + " is out of range.");
    T& element = _matrix[i][j];
    return element;
}

template <typename T>
T Linear::Matrix<T>::Sel(size_t i, size_t j) const {
    if (i >= _rows) throw OutOfRange("Row value " + std::to_string(i) + " is out of range.");
    if (j >= _columns) throw OutOfRange("Column value " + std::to_string(j) + " is out of range.");

    return _matrix[i][j];
}


template <typename T>
T Linear::Matrix<T>::Determinant() const {
    if (_columns != _rows)
        throw Linear::IncompatibleSize("Matrix is not a quadratic. Determinant cannot be found!");
    if (_columns == 0) return T(0);
    if (_columns == 1) return _matrix[0][0];

    T determinant = 0;
    for (size_t i = 0; i < _columns; ++i) {
        Matrix<T> minor(_rows - 1, _columns - 1);
        for (size_t j = 1; j < _rows; ++j) {
            size_t decrease = 0;
            for (size_t k = 0; k < _columns; ++k) {
                if (k == i) {
                    decrease = 1;
                    continue;
                }
                minor.At(j - 1, k - decrease) = _matrix[j][k];
            }

        }
        determinant += std::pow(T(-1), 2 + i) * _matrix[0][i] * minor.Determinant();
    }

    return determinant;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::Transpose() {
    std::vector<std::vector<T>> new_matrix(_columns, std::vector<T>(_rows, T()));
    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < _columns; ++j) {
            new_matrix[j][i] = _matrix[i][j];
        }
    }

    return Matrix<T>(new_matrix);
}

template <typename T>
std::pair<size_t, size_t> Linear::Matrix<T>::Size() const {
    return {_rows, _columns};
}


template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::operator+(const Matrix<T> &m) const {
    if (m.Size() != Size()) throw  Linear::IncompatibleSize("Cannot add matrixes with different sizes");

    Matrix<T> new_matrix(_rows, _columns);

    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < _columns; ++j) {
            new_matrix.At(i, j) = Sel(i, j) + m.Sel(i, j);
        }
    }

    return new_matrix;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::operator-(const Matrix<T> &m) const {
    return *this + T(-1) * m;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::operator*(T number) const {
    Matrix<T> new_matrix(_rows, _columns);

    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < _columns; ++j) {
            new_matrix.At(i, j) = Sel(i, j) * number;
        }
    }

    return new_matrix;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::dot(const Matrix<T> &m) const {
    auto size = m.Size();
    if (_columns != size.first) throw Linear::IncompatibleSize("Can't multiply matrix when columns of the first one (" +
                                                     std::to_string(_columns) +
                                                     ") do not equal to rows of the second (" +
                                                     std::to_string(size.first) + ")!");
    Matrix<T> new_matrix(_rows, size.second);
    for (size_t i = 0; i < _rows; ++i) {
        for (size_t j = 0; j < size.second; ++j) {
            for (size_t k = 0; k < _columns; ++k) {
                new_matrix.At(i, j) += Sel(i, k) * m.Sel(k, j);
            }
        }
    }

    return new_matrix;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::operator*(const Matrix<T> &m) const {
    if (Size() != m.Size()) throw Linear::IncompatibleSize("Sizes are incompatible, cannot do the "
                                                           "componentwise multiplication");

    Matrix<T> new_matrix(_rows, _columns);

    for (size_t i = 0; i < _rows; ++i){
        for (size_t j = 0; j < _columns; ++j) {
            new_matrix.At(i, j) = Sel(i, j) * m.Sel(i, j);
        }
    }

    return new_matrix;
}

template <typename T>
Linear::Diagonal<T>::Diagonal(const std::vector<std::vector<T>>& matrix) : Matrix<T>(matrix) {
    if (Matrix<T>::_rows != Matrix<T>::_columns) throw Linear::IncompatibleSize("Diagonal matrix "
                                                                               "must be quadratic");
    for (size_t i = 0; i < Matrix<T>::_rows; ++i){
        for (size_t j = 0; j < Matrix<T>::_columns; ++j) {
            if (i == j) continue;
            if (Matrix<T>::_matrix[i][j] != 0) throw Linear::WrongFormat("Diagonal matrix elements which does "
                                                                         "not belong to the diagonal must be zeros ");
        }
    }
}

template <typename T>
Linear::Diagonal<T>::Diagonal(const std::vector<T>& diagonal) : Matrix<T>(diagonal.size()) {
    for (size_t i = 0; i < Matrix<T>::_rows; ++i) {
        Matrix<T>::At(i, i) = diagonal[i];
    }
}

template <typename T>
T& Linear::Diagonal<T>::At(size_t i, size_t j) {
    if (i != j) throw Linear::OutOfRange("You can't change not diagonal elements in diagonal matrix");
    return Matrix<T>::At(i, j);
}

template <typename T>
T& Linear::Diagonal<T>::At(size_t i) {
    return Matrix<T>::At(i, i);
}

template <typename T>
T Linear::Diagonal<T>::Determinant() const {
    T determinant = 1;
    for (size_t i = 0; i < Matrix<T>::_rows; ++i) {
        determinant *= Matrix<T>::Sel(i, i);
    }

    return determinant;
}

template <typename T>
[[maybe_unused]] Linear::Diagonal<T>::Diagonal(int N) : Matrix<T>(N) {
    for (size_t i = 0; i < Matrix<T>::_rows; ++i) {
        Matrix<T>::At(i, i) = 1;
    }
}

template <typename T>
Linear::Unit<T>::Unit(int N) : Diagonal<T>(N) {}

template <typename T>
T Linear::Unit<T>::Determinant() const {
    return T(1);
}

template <typename T>
Linear::UpperTriangular<T>::UpperTriangular(const std::vector<std::vector<T>>& matrix) :Matrix<T>(matrix) {
    for (size_t i = 0; i < Matrix<T>::_rows; ++i) {
        for (size_t j = 0; j < i && j < Matrix<T>::_columns; ++j) {
            if (Matrix<T>::_matrix[i][j] != T(0)) throw Linear::OutOfRange("Upper triangular matrix cannot have "
                                                                           "non zero elements below the diagonal");
        }
    }
}

template <typename T>
T& Linear::UpperTriangular<T>::At(size_t i, size_t j) {
   if (i > j) throw Linear::OutOfRange("Cannot access the elements below the diagonal "
                                       "in the upper triangular matrix");
   return Matrix<T>::At(i, j);
}

template <typename T>
Linear::LowerTriangular<T>::LowerTriangular(const std::vector<std::vector<T>>& matrix) :Matrix<T>(matrix) {
    for (size_t i = 0; i < Matrix<T>::_rows; ++i) {
        for (size_t j = i + 1; j < Matrix<T>::_columns; ++j) {
            if (Matrix<T>::_matrix[i][j] != T(0)) throw Linear::OutOfRange("Lower triangular matrix cannot have "
                                                                           "non zero elements above the diagonal");
        }
    }
}

template <typename T>
T& Linear::LowerTriangular<T>::At(size_t i, size_t j) {
    if (i < j) throw Linear::OutOfRange("Cannot access the elements above the diagonal "
                                        "in the lower triangular matrix");
    return Matrix<T>::At(i, j);
}
