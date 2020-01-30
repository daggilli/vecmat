#ifndef __VECMAT_H__
#define __VECMAT_H__
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace VecMat {
  namespace {
    constexpr std::size_t GAUSS_ELIM_DET_THRESHOLD = 7;
  }

  template <typename T>
  class Matrix;

  template <typename T>
  class Vector {
   public:
    using iterator = T*;
    using const_iterator = T const*;
    using reverse_iterator = T*;
    using const_reverse_iterator = T const*;

    Vector() : vec(nullptr), sz(0), zero(0.0) {}
    Vector(const std::size_t n) : vec(new T[n]), sz(n), zero(0.0) {}
    Vector(const std::initializer_list<T> vlist) : vec(new T[vlist.size()]), sz(vlist.size()), zero(0.0) {
      std::size_t ix = 0;

      for (auto& el : vlist) {
        vec[ix++] = el;
      }
    }
    Vector(const Vector<T>& v) : vec(new T[v.sz]), sz(v.sz), zero(v.zero) { std::copy(v.vec, v.vec + sz, vec); }
    Vector(Vector<T>&& v) : vec(std::exchange(v.vec, nullptr)), sz(std::exchange(v.sz, 0)), zero(v.zero) {}
    ~Vector() { delete[] vec; }
    std::size_t size() const { return sz; }
    const T total() const { return std::accumulate(vec, vec + sz, zero); }
    const Vector<T> conj() const {
      Vector<T> v(this);
      return v;
    }
    const double norm() const {
      return std::sqrt(std::accumulate(vec, vec + sz, 0.0, [](const T& s, const T& v) { return s + v * v; }));
    }
    Vector<T>& operator=(const Vector<T>& v) {
      if (&v != this) {
        if (sz != v.sz) {
          delete[] vec;
          vec = new T[v.sz];
          sz = v.sz;
        }
        std::copy(v.vec, v.vec + sz, vec);
      }
      return *this;
    }
    Vector<T>& operator=(Vector<T>&& v) {
      if (&v != this) {
        vec = std::exchange(v.vec, nullptr);
        sz = std::exchange(v.sz, 0);
        zero = v.zero;
      }
      return *this;
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

      for (auto i = 0; i < lhs.size(); i++) v[i] = lhs[i] - rhs[i];

      return v;
    }
    Vector<T>& operator+=(const Vector<T>& rhs) {
      checkDimensions(*this, rhs);
      for (auto i = 0; i < vec.size(); i++) vec[i] += rhs.vec[i];
      return *this;
    }
    Vector<T>& operator-=(const Vector<T>& rhs) {
      checkDimensions(*this, rhs);
      for (auto i = 0; i < vec.size(); i++) vec[i] -= rhs.vec[i];
      return *this;
    }

    friend T operator*(const Vector<T>& lhs, const Vector<T>& rhs) {
      checkDimensions(lhs, rhs);
      T dot = 0;
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
      const_iterator it = v.cbegin();
      while (it != v.cend()) {
        os << sep << *it++;
        sep = " ";
      }
      os << "|";
      return os;
    }

    inline iterator begin() noexcept { return vec; }
    inline const_iterator cbegin() const noexcept { return vec; }
    inline iterator end() noexcept { return vec + sz; }
    inline const_iterator cend() const noexcept { return vec + sz; }
    inline reverse_iterator rbegin() noexcept { return vec + sz; }
    inline reverse_iterator crbegin() const noexcept { return vec + sz; }
    inline reverse_iterator rend() noexcept { return vec; }
    inline reverse_iterator crend() const noexcept { return vec; }

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
    T* vec;
    std::size_t sz;
    T zero;
  };

  template <>
  const inline Vector<std::complex<double>> Vector<std::complex<double>>::conj() const {
    Vector<std::complex<double>> v(sz);
    for (auto i = 0; i < sz; i++) {
      v[i] = std::conj(vec[i]);
    }
    return v;
  }

  template <>
  const inline double Vector<std::complex<double>>::norm() const {
    return std::sqrt(std::accumulate(
        vec, vec + sz, 0.0, [](const double s, const std::complex<double> t) { return s + std::abs(t * std::conj(t)); }));
  }

  template <typename T>
  const inline T angle(const Vector<T>& u, const Vector<T>& v) {
    return T(std::acos((u * v) / (u.norm() * v.norm())));
  }

  const inline std::complex<double> angle(const Vector<std::complex<double>>& u, const Vector<std::complex<double>>& v) {
    return std::acos(u * v.conj() / (u.norm() * v.norm()));
  }

  template <typename T>
  class Matrix {
    using iterator = T**;
    using const_iterator = const T**;
    using reverse_iterator = T**;
    using const_reverse_iterator = const T**;

   public:
    Matrix(const std::size_t n) : m(n), n(n), data(new T[n * n]), mat(new T*[n]) {
      for (auto i = 0; i < n; i++) mat[i] = data + i * n;
    }
    Matrix(const std::size_t m, const std::size_t n) : m(m), n(n), data(new T[m * n]), mat(new T*[m]) {
      std::cout << "rc ctor " << m << " x " << n << "\n";
      for (auto i = 0; i < m; i++) mat[i] = data + i * n;
    }
    Matrix(const std::pair<std::size_t, std::size_t>& sz)
        : m(sz.first), n(sz.second), data(new T[sz.first * sz.first]), mat(new T*[sz.first]) {
      std::cout << "rc sz ctor " << m << " x " << n << "\n";
      for (auto i = 0; i < m; i++) mat[i] = data + i * n;
    }
    Matrix(const std::initializer_list<std::initializer_list<T>> mlist)
        : m(mlist.size()), n(0), data(nullptr), mat(new T*[mlist.size()]) {
      std::size_t ix = 0;
      std::size_t jx;
      n = mlist.begin()->size();
      data = new T[m * n];
      std::cout << "il ctor " << m << " x " << n << "\n";
      for (auto& r : mlist) {
        jx = 0;
        for (auto& el : r) {
          data[(ix * n) + jx++] = el;
        }
        mat[ix++] = data + ix * n;
      }
    }
    ~Matrix() {
      delete[] mat;
      delete[] data;
    }
    const std::pair<std::size_t, std::size_t> size() const { return std::make_pair(m, n); }
    const std::size_t rows() const { return m; }
    const std::size_t cols() const { return n; }
    const bool isSquare() { return m == n; }
    const T trace() const {
      T tr = 0;
      for (std::size_t i = 0; i < std::min(m, n); i++) tr += mat[i][i];
      return tr;
    }
    const T determinant() const {
      checkSquare();
      x3 = 0;
      if (m >= GAUSS_ELIM_DET_THRESHOLD) return gaussDet(mat, m);
      auto d = rdet(mat, m);
      std::cout << "3x3: " << x3 << "\n";
      return d;
    }
    T* operator[](const std::size_t i) { return mat[i]; }
    const T* operator[](const std::size_t i) const { return mat[i]; }
    friend Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
      checkDimensions(lhs, rhs);
      auto lr = lhs.size().first;
      auto lc = lhs.size().second;
      auto rr = rhs.size().first;
      auto rc = rhs.size().second;

      std::cout << "lr " << lr << " lc " << lc << " rc " << rc << "\n";
      Matrix<T> ma(lr, rc);
      for (auto i = 0; i < lr; i++) {
        auto lv = lhs[i];
        for (auto j = 0; j < rr; j++) {
          auto rv = rhs[j];
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
        auto row = m[i];
        for (auto j = 0; j < mc; j++) {
          vec[i] += row[j] * v[j];
        }
      }
      return vec;
    }
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
      if (m.rows()) {
        std::string lf = "";
        const_iterator it = m.cbegin();
        while (it != m.cend()) {
          os << lf << "|";
          std::string sep = "";
          std::for_each(*it, *it + m.n, [&](const T& t) {
            os << sep << t;
            sep = " ";
          });
          os << "|";
          lf = '\n';
          it++;
        }
      } else
        os << "||";
      return os;
    }

    inline iterator begin() noexcept { return mat; }
    inline const_iterator cbegin() const noexcept { return const_cast<const double**>(mat); }
    inline iterator end() noexcept { return mat + m; }
    inline const_iterator cend() const noexcept { return const_cast<const double**>(mat + m); }
    inline reverse_iterator rbegin() noexcept { return mat + m; }
    inline reverse_iterator crbegin() const noexcept { return const_cast<const double**>(mat + m); }
    inline reverse_iterator rend() noexcept { return mat; }
    inline reverse_iterator crend() const noexcept { return const_cast<const double**>(mat); }

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
    static void checkEqualDimensions(const Matrix<T>& left, const Matrix<T>& right) {
      if (left.size().first != right.size().first || left.size().second != right.size().sexond) {
        throw(std::runtime_error("Matrix + matrix length mismatch"));
      }
    }
    const void checkSquare() const {
      if (m != n) {
        throw(std::runtime_error("Matrix not square"));
      }
    }
    const T rdet(T** const msub, const std::size_t sz) const {
      T det(0);
      if (sz < 1) {
        throw(std::runtime_error("Matrix empty"));
      }

      if (sz == 1)
        det = msub[0][0];
      else if (sz == 2)
        det = msub[0][0] * msub[1][1] - msub[1][0] * msub[0][1];
      else if (sz == 3) {
        x3++;
        det = (msub[0][0] * (msub[1][1] * msub[2][2] - msub[2][1] * msub[1][2])) -
              (msub[0][1] * (msub[1][0] * msub[2][2] - msub[1][2] * msub[2][0])) +
              (msub[0][2] * (msub[1][0] * msub[2][1] - msub[1][1] * msub[2][0]));
      } else {
        auto mind = sz - 1;
        T* d = new T[mind * mind];
        T** t = new T*[mind];
        for (auto i = 0; i < mind; i++) {
          t[i] = d + i * mind;
        }

        int_least16_t p = 1;
        for (auto i = 0; i < sz; i++) {
          getMinor(msub, t, 0, i, sz);
          det += rdet(t, sz - 1) * msub[0][i] * p;
          p = -p;
        }

        delete[] t;
        delete[] d;
      }

      return det;
    }

    const void getMinor(T** const m, T** const t, const std::size_t i, const std::size_t j, const std::size_t n) const {
      int ix, jx, ic = 0, jc;

      for (ix = 0; ix < n; ix++) {
        if (i == ix) continue;
        jc = 0;
        for (jx = 0; jx < n; jx++) {
          if (j == jx) continue;
          t[ic][jc] = m[ix][jx];
          jc++;
        }
        ic++;
      }
    }

    const T gaussDet(T** const msub, const std::size_t sz) const {
      T det(1);
      T* ad = new T[sz * sz];
      T** a = new T*[sz];
      std::copy(data, data + sz * sz, ad);
      for (auto i = 0; i < sz; i++) a[i] = ad + i * sz;

      uint parity = 0;

      for (auto i = 0; i < sz; i++) {
        for (auto k = i + 1; k < sz; k++) {
          if (std::abs(a[i][i]) < std::abs(a[k][i])) {
            parity ^= 1;
            for (auto j = 0; j < sz; j++) {
              std::swap(a[i][j], a[k][j]);
            }
          }
        }
      }

      for (auto i = 0; i < sz - 1; i++) {
        for (auto k = i + 1; k < sz; k++) {
          T t = a[k][i] / a[i][i];
          for (auto j = 0; j < sz; j++) {
            a[k][j] = a[k][j] - t * a[i][j];
          }
        }
      }

      for (auto i = 0; i < sz; i++) det *= a[i][i];

      return parity ? -det : det;
    }

    std::size_t m;
    std::size_t n;
    T* data;
    T** mat;
    mutable uint x3;
  };  // namespace VecMat

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

  template <typename T>
  const inline Matrix<T> conj(const Matrix<T>& m) {
    Matrix<T> ma(m.size());
    for (auto i = 0; i < m.size().second; i++) {
      for (auto j = 0; j < m.size().first; j++) {
        ma[i][j] = std::conj(m[i][j]);
      }
    }

    return ma;
  }

  template <typename T>
  const inline Matrix<T> hermitian(const Matrix<T>& m) {
    Matrix<T> ma(m.size().second, m.size().first);
    for (auto i = 0; i < m.size().second; i++) {
      for (auto j = 0; j < m.size().first; j++) {
        ma[j][i] = std::conj(m[i][j]);
      }
    }

    return ma;
  }
}  // namespace VecMat
#endif
