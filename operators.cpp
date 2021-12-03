#include "operators.h"

template <typename T>
std::ostream& operator<< (std::ostream& os, const Linear::Matrix<T>& matrix) {
    std::pair<size_t, size_t> size = matrix.Size();

    os << "[";
    for (size_t i = 0; i < size.first; ++i) {
        os << "[";
        for (size_t j = 0; j < size.second; ++j) {
            if (j != 0) {
                os << ", ";
            }
            os << matrix.Sel(i, j);
        }
        if (i != size.first - 1) {
            os << "],\n";
        } else {
            os << "]]";
        }
    }

    return os;
}

template <typename T>
Linear::Matrix<T> operator* (T number, const Linear::Matrix<T>& matrix) {
    auto size = matrix.Size();
    Linear::Matrix <T> new_matrix(size.first, size.second);

    for (int i = 0; i < size.first; ++i){
        for (int j = 0; j < size.second; ++j) {
            new_matrix.At(i, j) = matrix.Sel(i, j) * number;
        }
    }

    return new_matrix;
}
