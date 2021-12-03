#pragma once
#include <vector>

namespace Linear {
    template <typename T>
    class Matrix {
    public:
        explicit Matrix(const std::vector<std::vector<T>>&);
        explicit Matrix(const std::vector<T>&);
        Matrix(int, int);
        explicit Matrix(int);
        Matrix();

        virtual T& At(size_t i, size_t j);
        T Sel(size_t i, size_t j) const;
        [[nodiscard]] std::pair<size_t, size_t> Size() const;

        Matrix<T> Transpose();
        virtual T Determinant() const;

        Matrix<T> operator+ (const Matrix<T>& m) const;
        Matrix<T> operator- (const Matrix<T>& m) const;
        Matrix<T> operator* (T number) const;
        Matrix<T> dot(const Matrix<T>& m) const;
        Matrix<T> operator* (const Matrix<T>& m) const;
    protected:
        std::vector<std::vector<T>> _matrix;
        size_t _rows, _columns;
    };

    template <typename T>
    class Diagonal : public Matrix<T> {
    public:
        explicit Diagonal(const std::vector<std::vector<T>>&);
        explicit Diagonal(const std::vector<T>&);

        T& At(size_t i, size_t j) override;
        virtual T& At(size_t i);

        T Determinant() const override;

        [[maybe_unused]] explicit Diagonal(int);
    };

    template <typename T>
    class Unit : public Diagonal<T> {
    public:
        explicit Unit(int);

        T Determinant() const override;
    private:
        using Diagonal<T>::At;
    };

    template <typename T>
    class UpperTriangular : public Matrix<T> {
    public:
        explicit UpperTriangular(const std::vector<std::vector<T>>&);
        T& At(size_t i, size_t j) override;
    };

    template <typename T>
    class LowerTriangular : public Matrix<T> {
    public:
        explicit LowerTriangular(const std::vector<std::vector<T>>&);
        T& At(size_t i, size_t j) override;
    };
}
