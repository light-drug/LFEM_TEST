/*
 * Copyright (C) 2024 Xiamen University
 *
 * @Author: Liang Pan
 * @Date:   2024-12-02
 * @Last Modified by: Liang Pan
 * @Last Modified time: 2024-12-02
 */

#ifndef VECTOR_OVERLOAD_HPP
#define VECTOR_OVERLOAD_HPP

#include <stdexcept>
#include <vector>

template<typename T>
class VectorOperations {
public:
  // 实现数乘：std::vector<T> 和标量 Scalar
  template<typename Scalar>
  static std::vector<T> multiply(const std::vector<T>& vec, const Scalar& scalar) {
    const int& size = vec.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = vec[i] * scalar;  // 对于每个元素进行数乘
    }
    return result;
  };

  template<typename Scalar>
  static std::vector<T> multiply(const Scalar& scalar, const std::vector<T>& vec) {
    const int& size = vec.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = scalar * vec[i];  // 对于每个元素进行数乘
    }
    return result;
  };

  template<typename Scalar>
  static std::vector<T> divide(const std::vector<T>& vec, const Scalar& scalar) {
    const int& size = vec.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = vec[i] / scalar;  // 对于每个元素进行数乘
    }
    return result;
  };

  // 实现相加
  static std::vector<T> add(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors must be of the same size to add.");
    }
    const int& size = vec1.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = vec1[i] + vec2[i];  // 对应元素相加
    }
    return result;
  };

  // 实现相减
  static std::vector<T> subtract(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors must be of the same size to subtract.");
    }
    const int& size = vec1.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = vec1[i] - vec2[i];  // 对应元素相减
    }
    return result;
  };
};

// 重载 * 运算符，使得 std::vector<T> 和标量相乘
template<typename T, typename Scalar>
std::vector<T> operator*(const std::vector<T>& vec, const Scalar& scalar) {
  return VectorOperations<T>::multiply(vec, scalar);
};

// 重载 * 运算符，使得标量和 std::vector<T> 相乘（不支持交换顺序）
template<typename T, typename Scalar>
std::vector<T> operator*(const Scalar& scalar, const std::vector<T>& vec) {
  return VectorOperations<T>::multiply(scalar, vec);
};

// 重载 + 运算符，使得两个相同类型的 std::vector 相加
template<typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2) {
  return VectorOperations<T>::add(vec1, vec2);
}

template<typename T>
std::vector<T>& operator+=(std::vector<T>& vec1, const std::vector<T>& vec2) {
    vec1 = VectorOperations<T>::add(vec1, vec2); // 更新 vec1
    return vec1; // 返回更新后的 vec1 引用
}

// 重载 - 运算符，使得两个相同类型的 std::vector 相减
template<typename T>
std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2) {
  return VectorOperations<T>::subtract(vec1, vec2);
}

template<typename T, typename Scalar>
std::vector<T> operator/(const std::vector<T>& vec, const Scalar& scalar) {
  return VectorOperations<T>::divide(vec, scalar);
};

#endif  //  VECTOR_OVERLOAD 
