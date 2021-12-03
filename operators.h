#pragma once

#include <iostream>
#include "linear.h"

template <typename T>
std::ostream& operator<< (std::ostream& os, const Linear::Matrix<T>& matrix);

template <typename T>
Linear::Matrix<T> operator* (T number, const Linear::Matrix<T>& matrix);