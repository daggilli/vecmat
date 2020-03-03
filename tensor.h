#ifndef __TENSOR_H__
#define __TENSOR_H__
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace Tensor {
  namespace {
    constexpr std::size_t GAUSS_ELIM_DET_THRESHOLD = 7;
  }

  template <typename T>
  class Matrix;

  template <typename T>
  using VectorGenerator = std::function<T(const std::size_t index)>;

  template <typename T>
  class VectorGeneratorFunctor {
   public:
    virtual const T operator()(const std::size_t index) const = 0;
  };

  template <typename T>
  class Vector {
   public:
    using iterator = T*;
    using const_iterator = T const*;
    using reverse_iterator = T*;
    using const_reverse_iterator = T const*;

    Vector() : vec(nullptr), sz(0), zero(0.0) {}
    Vector(const std::size_t n) : vec(new T[n]), sz(n), zero(0.0) {}
    Vector(const std::size_t n, VectorGenerator<T> gen) : vec(new T[n]()), sz(n), zero(0.0) {
      for (auto i = 0; i < n; i++) vec[i] = gen(i);
    }

    Vector(const std::initializer_list<T> vlist) : vec(new T[vlist.size()]), sz(vlist.size()), zero(0.0) {
      std::size_t ix = 0;

      for (auto& el : vlist) {
        vec[ix++] = el;
      }
    }
    Vector(const Vector<T>& v) : vec(new T[v.sz]), sz(v.sz), zero(v.zero) { std::copy(v.vec, v.vec + sz, vec); }
    Vector(Vector<T>&& v) : vec(std::exchange(v.vec, nullptr)), sz(std::exchange(v.sz, 0)), zero(v.zero) {}
    virtual ~Vector() { delete[] vec; }
    T& operator[](const std::size_t i) { return vec[i]; }
    const T& operator[](const std::size_t i) const { return vec[i]; }
    T& operator()(const std::size_t i) { return vec[i]; }
    const T& operator()(const std::size_t i) const { return vec[i]; }
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
      for (auto i = 0; i < sz; i++) vec[i] += rhs.vec[i];
      return *this;
    }
    Vector<T>& operator-=(const Vector<T>& rhs) {
      checkDimensions(*this, rhs);
      for (auto i = 0; i < sz; i++) vec[i] -= rhs.vec[i];
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
    std::size_t size() const { return sz; }
    const T total() const { return std::accumulate(vec, vec + sz, zero); }
    const Vector<T> conj() const {
      Vector<T> v(this);
      return v;
    }
    const double norm() const {
      return std::sqrt(std::accumulate(vec, vec + sz, 0.0, [](const T& s, const T& v) { return s + v * v; }));
    }
    const T mean() const { return std::accumulate(vec, vec + sz, zero) / static_cast<double>(sz); }
    const double variance() const {
      auto avg = mean();
      return std::accumulate(vec, vec + sz, 0.0,
                             [&avg](const double s, const T& v) {
                               auto diff = std::abs(avg - v);
                               return s + diff * diff;
                             }) /
             static_cast<double>(sz - 1);
    }
    const double stddev() const { return std::sqrt(variance()); }
    const T rms() const {
      return std::sqrt(std::accumulate(vec, vec + sz, zero, [](const T& s, const T& v) { return s + v * v; }) /
                       static_cast<double>(sz));
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
      if (!checkDimensionsNoThrow(left, right)) {
        throw(std::runtime_error("Vector x vector length mismatch"));
      }
    }
    static void checkDimensions(const Vector<T>& left, const Matrix<T>& right) {
      if (!checkDimensionsNoThrow(left, right)) {
        throw(std::runtime_error("Vector x matrix length mismatch"));
      }
    }
    static bool checkDimensionsNoThrow(const Vector<T>& left, const Vector<T>& right) { return left.size() == right.size(); }
    static bool checkDimensionsnoThrow(const Vector<T>& left, const Matrix<T>& right) {
      return left.size() == right.size().first;
    }
    T* vec;
    std::size_t sz;
    T zero;
  };  // class Vector

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

  template <typename T>
  const inline std::complex<T> angle(const Vector<std::complex<T>>& u, const Vector<std::complex<T>>& v) {
    return std::acos(u * v.conj() / (u.norm() * v.norm()));
  }

  template <typename T>
  using MatrixGenerator = std::function<T(const std::size_t, const std::size_t)>;

  template <typename T>
  class MatrixGeneratorFunctor {
   public:
    virtual const T operator()(const std::size_t i, const std::size_t j) const = 0;
  };

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
    Matrix(const std::size_t n, MatrixGenerator<T> gen) : m(n), n(n), data(new T[n * n]), mat(new T*[n]) {
      std::size_t ix = 0;
      for (auto i = 0; i < n; i++) {
        mat[i] = data + i * n;
        for (auto j = 0; j < n; j++) data[ix++] = gen(i, j);
      }
    }
    Matrix(const std::size_t m, const std::size_t n) : m(m), n(n), data(new T[m * n]), mat(new T*[m]) {
      for (auto i = 0; i < m; i++) mat[i] = data + i * n;
    }
    Matrix(const std::size_t m, const std::size_t n, MatrixGenerator<T> gen) : m(m), n(n), data(new T[m * n]), mat(new T*[m]) {
      std::size_t ix = 0;
      for (auto i = 0; i < m; i++) {
        mat[i] = data + i * n;
        for (auto j = 0; j < n; j++) data[ix++] = gen(i, j);
      }
    }
    Matrix(const std::pair<std::size_t, std::size_t>& sz)
        : m(sz.first), n(sz.second), data(new T[sz.first * sz.first]), mat(new T*[sz.first]) {
      for (auto i = 0; i < m; i++) mat[i] = data + i * n;
    }
    Matrix(const std::pair<std::size_t, std::size_t>& sz, MatrixGenerator<T> gen)
        : m(sz.first), n(sz.second), data(new T[sz.first * sz.first]), mat(new T*[sz.first]) {
      std::size_t ix = 0;
      for (auto i = 0; i < m; i++) {
        mat[i] = data + i * n;
        for (auto j = 0; j < n; j++) data[ix++] = gen(i, j);
      }
    }
    Matrix(const std::initializer_list<std::initializer_list<T>> mlist)
        : m(mlist.size()), n(0), data(nullptr), mat(new T*[mlist.size()]) {
      std::size_t ix = 0;
      std::size_t jx;
      n = mlist.begin()->size();
      data = new T[m * n];

      for (auto& r : mlist) {
        jx = 0;
        for (auto& el : r) {
          data[(ix * n) + jx++] = el;
        }
        mat[ix++] = data + ix * n;
      }
    }
    virtual ~Matrix() {
      delete[] mat;
      delete[] data;
    }
    const std::pair<std::size_t, std::size_t> size() const { return std::make_pair(m, n); }
    T* operator[](const std::size_t i) { return mat[i]; }
    const T* operator[](const std::size_t i) const { return mat[i]; }
    T& operator()(const std::size_t i, const std::size_t j) { return data[i * n + j]; }
    const T& operator()(const std::size_t i, const std::size_t j) const { return data[i * n + j]; }
    friend Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
      checkEqualDimensions(lhs, rhs);
      Matrix<T> ma(lhs.size());

      for (auto i = 0; i < lhs.size().first; i++) {
        for (auto j = 0; j < lhs.size().second; j++) {
          ma[i][j] = lhs[i][j] + rhs[i][j];
        }
      }
      return ma;
    }
    friend Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
      checkDimensions(lhs, rhs);
      auto lr = lhs.size().first;
      auto lc = lhs.size().second;
      auto rr = rhs.size().first;
      auto rc = rhs.size().second;

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
    friend Matrix<T> operator*(const Matrix<T>& lhs, const T scal) {
      Matrix<T> ma(lhs.size());

      for (auto i = 0; i < lhs.size().first; i++) {
        for (auto j = 0; j < lhs.size().second; j++) {
          ma[i][j] = lhs[i][j] * scal;
        }
      }

      return ma;
    }
    friend Matrix<T> operator*(const T scal, const Matrix<T>& lhs) { return lhs * scal; }
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
      if (m >= GAUSS_ELIM_DET_THRESHOLD) return gaussDet();
      return rdet(mat, m);
    }
    const Matrix<T> inverse() const {
      if (m < 4) return invert();
      return gaussInvert();
    }

    inline iterator begin() noexcept { return mat; }
    inline const_iterator cbegin() const noexcept { return const_cast<const T**>(mat); }
    inline iterator end() noexcept { return mat + m; }
    inline const_iterator cend() const noexcept { return const_cast<const T**>(mat + m); }
    inline reverse_iterator rbegin() noexcept { return mat + m; }
    inline reverse_iterator crbegin() const noexcept { return const_cast<const T**>(mat + m); }
    inline reverse_iterator rend() noexcept { return mat; }
    inline reverse_iterator crend() const noexcept { return const_cast<const T**>(mat); }

   private:
    static void checkDimensions(const Matrix<T>& left, const Matrix<T>& right) {
      if (!checkDimensionsNoThrow(left, right)) {
        throw(std::runtime_error("Matrix x matrix length mismatch"));
      }
    }
    static void checkDimensions(const Matrix<T>& m, const Vector<T>& v) {
      if (!checkDimensionsNoThrow(m, v)) {
        throw(std::runtime_error("Matrix x vector length mismatch"));
      }
    }
    static void checkEqualDimensions(const Matrix<T>& left, const Matrix<T>& right) {
      if (!checkEqualDimensionsNoThrow(left, right)) {
        throw(std::runtime_error("Matrix + matrix length mismatch"));
      }
    }
    const void checkSquare() const {
      if (!checkSquareNoThrow()) {
        throw(std::runtime_error("Matrix not square"));
      }
    }
    static bool checkDimensionsNoThrow(const Matrix<T>& left, const Matrix<T>& right) { return left.cols() == right.rows(); }
    static bool checkDimensionsNoThrow(const Matrix<T>& m, const Vector<T>& v) { return m.cols() == v.size(); }
    static bool checkEqualDimensionsNoThrow(const Matrix<T>& left, const Matrix<T>& right) {
      return left.rows() == right.rows() && left.cols() == right.cols();
    }
    const bool checkSquareNoThrow() const { return m == n; }
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

    const T gaussDet() const {
      T det(1);
      T* ad = new T[m * m];
      T** a = new T*[m];
      std::copy(data, data + m * m, ad);
      for (auto i = 0; i < m; i++) a[i] = ad + i * m;

      uint parity = 0;

      for (auto i = 0; i < m; i++) {
        for (auto k = i + 1; k < m; k++) {
          if (std::abs(a[i][i]) < std::abs(a[k][i])) {
            parity ^= 1;
            for (auto j = 0; j < m; j++) {
              std::swap(a[i][j], a[k][j]);
            }
          }
        }
      }

      for (auto i = 0; i < m - 1; i++) {
        for (auto k = i + 1; k < m; k++) {
          T t = a[k][i] / a[i][i];
          for (auto j = 0; j < m; j++) {
            a[k][j] = a[k][j] - t * a[i][j];
          }
        }
      }

      for (auto i = 0; i < m; i++) det *= a[i][i];

      delete[] a;
      delete[] ad;
      return parity ? -det : det;
    }

    const Matrix<T> invert() const {
      Matrix<T> ma(size());

      if (m == 1) {
        ma.data[0] = 1 / data[0];
      } else if (m == 2) {
        T idet = 1 / (data[0] * data[3] - data[1] * data[2]);
        ma.data[0] = data[3] * idet;
        ma.data[1] = -data[1] * idet;
        ma.data[2] = -data[2] * idet;
        ma.data[3] = data[0] * idet;
      } else if (m == 3) {
        T idet = 1 / (data[0] * (data[4] * data[8] - data[5] * data[7]) + data[1] * (data[5] * data[6] - data[3] * data[8]) +
                      data[2] * (data[3] * data[7] - data[4] * data[6]));
        ma.data[0] = idet * (data[4] * data[8] - data[5] * data[7]);
        ma.data[1] = idet * (data[2] * data[7] - data[1] * data[8]);
        ma.data[2] = idet * (data[1] * data[5] - data[2] * data[4]);
        ma.data[3] = idet * (data[5] * data[6] - data[3] * data[8]);
        ma.data[4] = idet * (data[0] * data[8] - data[2] * data[6]);
        ma.data[5] = idet * (data[2] * data[3] - data[0] * data[5]);
        ma.data[6] = idet * (data[3] * data[7] - data[4] * data[6]);
        ma.data[7] = idet * (data[1] * data[6] - data[0] * data[7]);
        ma.data[8] = idet * (data[0] * data[4] - data[1] * data[3]);
      }

      return ma;
    }

    const Matrix<T> gaussInvert() const {
      Matrix<T> ma(size());

      std::size_t wd = 2 * m;
      T* ad = new T[m * wd]();
      T** a = new T*[m];
      for (auto i = 0; i < m; i++) {
        a[i] = ad + i * wd;
        std::copy(data + i * m, data + i * m + m, ad + i * wd);
        ad[i * wd + m + i] = T(1);
      }

      for (auto i = m - 1; i > 0; i--) {
        if (a[i - 1][0] < a[i][0]) {
          std::swap(a[i], a[i - 1]);
        }
      }

      for (auto i = 0; i < m; i++) {
        for (auto j = 0; j < m; j++) {
          if (i != j) {
            T t = a[j][i] / a[i][i];
            for (auto k = 0; k < wd; k++) {
              a[j][k] -= a[i][k] * t;
            }
          }
        }
      }

      for (auto i = 0; i < m; i++) {
        T t = a[i][i];
        for (auto j = m; j < wd; j++) {
          ma.mat[i][j - m] = a[i][j] / t;
        }
      }

      return ma;
    }

#ifndef NODEBUG
    void dump(T* d, const std::size_t r, const std::size_t c) const {
      std::string sep, lf;
      for (auto i = 0; i < r; i++) {
        sep = "";
        std::cout << lf << "|";
        for (auto j = 0; j < c; j++) {
          std::cout << sep << d[i * c + j];
          sep = " ";
        }
        std::cout << "|";
        lf = "\n";
      }
      std::cout << '\n';
    }

    void dump(T** d, const std::size_t r, const std::size_t c) const {
      std::string sep, lf;
      for (auto i = 0; i < r; i++) {
        sep = "";
        std::cout << lf << "|";
        for (auto j = 0; j < c; j++) {
          std::cout << sep << d[i][j];
          sep = " ";
        }
        std::cout << "|";
        lf = "\n";
      }
      std::cout << '\n';
    }
#endif

    std::size_t m;
    std::size_t n;
    T* data;
    T** mat;
    mutable uint x3;
  };  // class Matrix

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
}  // namespace Tensor
#endif
