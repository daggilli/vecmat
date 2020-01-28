#ifndef __VECMAT_H__
#define __VECMAT_H__
#include <cmath>
#include <complex>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace VecMat {
  template <typename T>
  class Matrix;

  template <typename T>
  class Vector {
   public:
    using vectorType = typename std::vector<T>;
    using iterator = typename vectorType::iterator;
    using const_iterator = typename vectorType::const_iterator;
    using reverse_iterator = typename vectorType::reverse_iterator;
    using const_reverse_iterator = typename vectorType::const_reverse_iterator;

    Vector() : zero(0.0) {}
    Vector(const std::size_t n) : vec(n), zero(0.0) {}
    Vector(const std::initializer_list<T> vlist) : vec(vlist), zero(0.0) {}

    void reserve(const std::size_t n) { vec.reserve(n); }
    void resize(const std::size_t n) { vec.resize(n); }
    std::size_t size() const { return vec.size(); }
    const T total() const { return std::accumulate(vec.cbegin(), vec.cend(), zero); }
    const Vector<T> conj() const {
      Vector<T> v = this;
      return v;
    }
    const double norm() const {
      return std::sqrt(std::accumulate(vec.cbegin(), vec.cend(), 0.0, [](const T& s, const T& v) { return s + v * v; }));
    }
    T& operator[](const std::size_t i) { return vec[i]; }
    const T& operator[](const std::size_t i) const { return vec[i]; }
    friend Vector<T> operator+(const Vector<T>& lhs, const Vector<T>& rhs) {
      checkDimensions(lhs, rhs);
      Vector<T> v(lhs.size());

      for (auto i = 0; i < lhs.size(); i++) v[i] = lhs[i] + rhs[i];

      return v;
    }
    friend Vector<T> operator-(const Vector<T>& lhs, const Vector<T>& rhs) {
      checkDimensions(lhs, rhs);
      Vector<T> v(lhs.size());

      for (auto i = 0; i < lhs.size(); i++) v[i] = lhs[i] + rhs[i];

      return v;
    }
    Vector<T>& operator+=(const Vector<T>& rhs) {
      checkDimensions(*this, rhs);
      for (auto i = 0; i < vec.size(); i++) vec[i] += rhs.vec[i];
      return *this;
    }
    Vector<T>& operator-=(const Vector<T>& rhs) {
      checkDimensions(*this, rhs);
      for (auto i = 0; i < vec.size(); i++) vec[i] += rhs.vec[i];
      return *this;
    }
    friend T operator*(const Vector<T>& lhs, const Vector<T>& rhs) {
      checkDimensions(lhs, rhs);
      T dot;
      for (auto i = 0; i < lhs.size(); i++) {
        dot += lhs[i] * rhs[i];
      }

      return dot;
    }
    friend Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m) {
      checkDimensions(v, m);
      std::size_t mr = m.size().first;
      std::size_t mc = m.size().second;

      Vector<T> nv(mc);

      for (auto i = 0; i < mr; i++) {
        auto row = m[i];
        for (auto j = 0; j < mc; j++) {
          nv[j] += v[i] * row[j];
        }
      }
      return nv;
    }
    friend Vector<T> operator*(const Vector<T>& v, const T& x) {
      std::size_t vs = v.size();
      Vector<T> nv(vs);

      for (auto i = 0; i < vs; i++) {
        nv[i] += v[i] * x;
      }

      return nv;
    }
    friend Vector<T> operator*(const T& x, const Vector<T>& v) { return v * x; }
    friend std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
      os << "|";
      std::string sep = "";
      for (auto& e : v.vec) {
        os << sep << e;
        sep = " ";
      }
      os << "|";
      return os;
    }

    inline iterator begin() noexcept { return vec.begin(); }
    inline const_iterator cbegin() const noexcept { return vec.cbegin(); }
    inline iterator end() noexcept { return vec.end(); }
    inline const_iterator cend() const noexcept { return vec.cend(); }
    inline reverse_iterator rbegin() noexcept { return vec.rbegin(); }
    inline reverse_iterator crbegin() const noexcept { return vec.crbegin(); }
    inline reverse_iterator rend() noexcept { return vec.rend(); }
    inline reverse_iterator crend() const noexcept { return vec.crend(); }

   private:
    static void checkDimensions(const Vector<T>& left, const Vector<T>& right) {
      if (left.size() != right.size()) {
        throw(std::runtime_error("Vector x vector length mismatch"));
      }
    }
    static void checkDimensions(const Vector<T>& left, const Matrix<T>& right) {
      if (left.size() != right.size().first) {
        throw(std::runtime_error("Vector x matrix length mismatch"));
      }
    }
    std::vector<T> vec;
    T zero;
  };

  template <>
  const inline Vector<std::complex<double>> Vector<std::complex<double>>::conj() const {
    Vector<std::complex<double>> v(vec.size());
    for (auto i = 0; i < vec.size(); i++) {
      v[i] = std::conj(vec[i]);
    }
    return v;
  }

  template <>
  const inline double Vector<std::complex<double>>::norm() const {
    return std::sqrt(std::accumulate(vec.cbegin(), vec.cend(), 0.0, [](const double s, const std::complex<double> t) {
      return s + std::abs(t * std::conj(t));
    }));
  }

  template <typename T>
  const inline double angle(const Vector<T>& u, const Vector<T>& v) {
    return std::acos((u * v) / (u.norm() * v.norm()));
  }

  const inline std::complex<double> angle(const Vector<std::complex<double>>& u, const Vector<std::complex<double>>& v) {
    return std::acos(u * v.conj() / (u.norm() * v.norm()));
  }

  template <typename T>
  class Matrix {
    using matrixType = typename std::vector<Vector<T>>;
    using iterator = typename matrixType::iterator;
    using const_iterator = typename matrixType::const_iterator;
    using reverse_iterator = typename matrixType::reverse_iterator;
    using const_reverse_iterator = typename matrixType::const_reverse_iterator;

   public:
    Matrix(const std::size_t n) : m(n), n(n), mat(n) {
      for (auto i = 0; i < n; i++) mat[i].resize(n);
    }
    Matrix(const std::size_t m, const std::size_t n) : m(m), n(n), mat(m) {
      std::cout << "rc ctor " << m << " x " << n << "\n";
      mat.resize(m);
      for (auto i = 0; i < m; i++) mat[i].resize(n);
    }
    Matrix(const std::initializer_list<std::initializer_list<T>> mlist) : m(mlist.size()), n(0), mat(mlist.size()) {
      auto mbeg = mat.begin();
      n = mlist.begin()->size();
      std::cout << "il ctor " << m << " x " << n << "\n";
      for (auto& v : mlist) {
        *mbeg++ = v;
      }
    }
    const std::pair<std::size_t, std::size_t> size() const { return std::make_pair(m, n); }
    const std::size_t rows() const { return m; }
    const std::size_t cols() const { return n; }
    const T trace() const {
      T tr = 0;
      for (std::size_t i = 0; i < mat.size(); i++) tr += mat[i][i];
      return tr;
    }
    Vector<T>& operator[](const std::size_t i) { return mat[i]; }
    const Vector<T>& operator[](const std::size_t i) const { return mat[i]; }
    friend Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
      checkDimensions(lhs, rhs);
      auto lr = lhs.size().first;
      auto lc = lhs.size().second;
      auto rr = rhs.size().first;
      auto rc = rhs.size().second;

      std::cout << "lr " << lr << " lc " << lc << " rc " << rc << "\n";
      Matrix<T> ma(lr, rc);
      std::cout << "++---+ " << ma[0].size() << "\n";
      for (auto i = 0; i < lr; i++) {
        auto& lv = lhs[i];
        for (auto j = 0; j < rr; j++) {
          auto& rv = rhs[j];
          for (auto k = 0; k < rc; k++) {
            ma[i][k] += lv[j] * rv[k];
          }
        }
      }

      return ma;
    }
    friend Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v) {
      std::cout << "op * m x v " << m.size().first << " " << m.size().second << " | " << v.size() << '\n';
      checkDimensions(m, v);
      std::size_t mr = m.size().first;
      std::size_t mc = v.size();

      Vector<T> vec(mr);
      for (auto i = 0; i < mr; i++) {
        auto& row = m[i];
        for (auto j = 0; j < mc; j++) {
          vec[i] += row[j] * v[j];
        }
      }
      return vec;
    }
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
      if (m.mat.size()) {
        std::string lf = "";
        for (auto& r : m.mat) {
          os << lf << r;
          lf = '\n';
        }
      } else
        os << "||";
      return os;
    }

    inline iterator begin() noexcept { return mat.begin(); }
    inline const_iterator cbegin() const noexcept { return mat.cbegin(); }
    inline iterator end() noexcept { return mat.end(); }
    inline const_iterator cend() const noexcept { return mat.cend(); }
    inline reverse_iterator rbegin() noexcept { return mat.rbegin(); }
    inline reverse_iterator crbegin() const noexcept { return mat.crbegin(); }
    inline reverse_iterator rend() noexcept { return mat.rend(); }
    inline reverse_iterator crend() const noexcept { return mat.crend(); }

   private:
    static void checkDimensions(const Matrix<T>& left, const Matrix<T>& right) {
      if (left.size().second != right.size().first) {
        throw(std::runtime_error("Matrix x matrix length mismatch"));
      }
    }
    static void checkDimensions(const Matrix<T>& m, const Vector<T>& v) {
      if (m.size().second != v.size()) {
        throw(std::runtime_error("Matrix x vector length mismatch"));
      }
    }
    std::size_t m;
    std::size_t n;
    std::vector<Vector<T>> mat;
  };

  template <typename T>
  const inline Matrix<T> transpose(const Matrix<T>& m) {
    Matrix<T> t(m.size().second, m.size().first);
    for (auto i = 0; i < m.size().second; i++) {
      for (auto j = 0; j < m.size().first; j++) {
        t[i][j] = m[j][i];
      }
    }

    return t;
  }
}  // namespace VecMat
#endif
