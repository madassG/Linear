#include <iostream>
#include "linear.cpp"

int main() {
    Linear::Matrix<int> matrix({{6, 7, 2}, {4, 6, 8}, {90, 13, 3}});
    Linear::Matrix<int> matrixx({{1, 2, 1}, {1, 1, 1}, {1, 1, 1}});
    Linear::LowerTriangular<int> m({{1, 0, 0}, {0, 5, 0}, {0, 0, 9}, {10, 11, 12}});
    std::cout << m << '\n';
    std::cout << matrix - matrixx << '\n';
    std::cout << Linear::Matrix<int>({{1, 4}, {2, 5}, {3, 6}}).dot(Linear::Matrix<int>({{7, 8, 9}, {10, 11, 12}})) << '\n';

    std::cout << matrix * matrixx << '\n';
    Linear::Unit<double> m1(10);
    std::cout << m1.Determinant() << '\n';

    Linear::Matrix<int> matrix2(5, 5);
    Linear::Matrix<int> matrix3(10);
    Linear::Matrix<int> matrix4;
    std::cout << matrix2.Size().first << ' ' << matrix2.Size().second << std::endl;
    matrix.Transpose();
    std::cout << matrix4 << '\n';
    std::cout << matrix4.Size().first << ' ' << matrix4.Size().second << std::endl;
    return 0;
}
