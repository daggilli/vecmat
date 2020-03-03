#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "tensor.h"

using namespace std::complex_literals;

double gen(const std::size_t i);

template <typename T>
class DerivedFunctor : public Tensor::MatrixGeneratorFunctor<T> {
 public:
  const T operator()(const std::size_t i, const std::size_t j) const {
    return std::sqrt(std::pow(static_cast<double>(i + 3), (static_cast<double>(j) - 1)));
  }
};

class DoubleFunctor : public Tensor::MatrixGeneratorFunctor<double> {
 public:
  DoubleFunctor() : count(0) {}
  const double operator()(const std::size_t i, const std::size_t j) const {
    count++;
    return count;
  }

 private:
  mutable std::size_t count;
};

template <typename T>
class OrderedPair : private std::pair<T, T> {
 public:
  OrderedPair() {}
  OrderedPair(std::pair<T, T> p) : std::pair<T, T>(p) {}
  using std::pair<T, T>::first;
  using std::pair<T, T>::second;
  friend std::ostream& operator<<(std::ostream& os, const OrderedPair<T>& v) {
    os << "[" << v.first << ", " << v.second << "]";
    return os;
  }
};

int main() {
  /*
  Tensor::Vector<float> fv(3);

  Tensor::Matrix<double> dm(3);
  dm[2][0] = 1.8;

  std::cout << dm[0][1] << '\n';

  Tensor::Matrix<double> dom(3, 2);

  Tensor::Vector<double> dva = {3.0, 4, 7, 1.9}, dvb = {6.7, 1.1, 5.3, 6.1};
  std::cout << dva.size() << '\n';
  Tensor::Vector<double> dvc = dva + dvb;

  for (auto e : dvc) {
    std::cout << e << " ";
  }
  std::cout << '\n';

  Tensor::Matrix mil = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  std::cout << dvc << '\n';

  std::cout << mil << '\n';

  dvc += dvc;
  std::cout << dvc << '\n';

  Tensor::Matrix mila = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Tensor::Matrix milb = {{10, 11, 12}, {13, 14, 15}, {16, 17, 18}};
  Tensor::Matrix milc = mila * milb;
  std::cout << milc << '\n';

  Tensor::Matrix a = {{1, 2, 3}, {4, 5, 6}};
  Tensor::Matrix b = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
  Tensor::Matrix c = a * b;
  std::cout << a << "\n\n";
  std::cout << b << "\n\n";
  std::cout << c << "\n"; */

  /*
    Tensor::Vector<double> v = {1, 2, 3, 4};
    Tensor::Vector<double> w = {1, 2, 3};
    Tensor::Matrix<double> n = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
    Tensor::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Tensor::Vector<double> r = n * v;

    std::cout << r << "\n";

    r = w * n;

    std::cout << r << "\n";

    std::cout << w * w << '\n';

    std::cout << s * n << '\n';

    std::cout << Tensor::transpose(n) << '\n';

    auto x = s * n;

    std::cout << x << '\n';

    auto y = std::move(x);

    std::cout << y << '\n';
    std::cout << "mv\n";
    std::cout << x << '\n';

    Tensor::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};

    std::cout << cv << '\n';

    std::cout << cv.total() << '\n';

    std::cout << w.norm() << '\n';

    std::cout << cv.norm() << '\n';

    std::cout << w * 3.3 << '\n';
    std::cout << 6.2 * w << '\n';

    v = {1, 2, 3};
    w = {4, 5, 6};

    std::cout << Tensor::angle(v, w) << '\n';

    cv = {1.0 + 1.0i, 2.0 - 1.0i};
    Tensor::Vector<std::complex<double>> cw = {2.0 - 3.0i, 4.0 + 1.0i};

    std::cout << cv.norm() << " " << cw.norm() << '\n';

    std::cout << Tensor::angle(cv, cw) << '\n';

    Tensor::Matrix<std::complex<double>> cm = {{1.0 + 1.0i, 2.0 - 1.0i}, {2.0 - 3.0i, 4.0 + 1.0i}};

    std::cout << cm.trace() << '\n';

    */

#if 0
  Tensor::Vector<double> v = {1, 2, 3, 4};
  Tensor::Vector<double> w = {1, 2, 3};
  Tensor::Matrix<double> n = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
  Tensor::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  // Tensor::Vector<double> r = n * v;
  auto r = w * n;

  std::cout << r << "\n---\n";
  std::cout << s[1][2] << '\n';

  // s[1][1] = 17;

  // std::cout << s[1][1] << '\n';

  std::cout << s << '\n';

  std::cout << s * n << '\n';

  auto sxn = s * n;
  std::cout << sxn.trace() << '\n';
  std::cout << w * w << '\n';

  std::cout << n * v << '\n';

  Tensor::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};

  std::cout << cv << '\n';

  std::cout << cv.total() << '\n';

  w = v;

  std::cout << w << '\n';

  Tensor::Vector<double> vv(std::move(w));

  std::cout << vv << '\n';

  w = std::move(vv);

  std::cout << w << '\n';

  std::cout << w.norm() << '\n';

  cv = {1.0 + 1.0i, 2.0 - 1.0i};
  Tensor::Vector<std::complex<double>> cw = {2.0 - 3.0i, 4.0 + 1.0i};

  std::cout << cv.norm() << " " << cw.norm() << '\n';

  std::cout << Tensor::angle(cv, cw) << '\n';

  std::cout << std::boolalpha << s.isSquare() << '\n';
  std::cout << std::boolalpha << sxn.isSquare() << '\n';

  s[1][1] = 17;

  Tensor::Matrix<double> g = {{3, 4, 22}, {14, 12, 30}, {26, 24, 42}};
  std::cout << g.determinant() << '\n';

  Tensor::Matrix<double> f = {{3, 4, 22, 4, 54, 4, 153},          {14, 12, 30, 12, 62, 18, 165},
                              {26, 24, 42, 24, 111, 36, 183},     {42, 40, 58, 60, 135, 60, 207},
                              {62, 60, 117, 90, 165, 90, 237},    {86, 126, 153, 126, 201, 126, 273},
                              {171, 168, 195, 168, 243, 168, 315}};

  std::cout << f.determinant() << '\n';

  Tensor::Matrix<int64_t> h = {
      {3, 4, 22, 4, 54}, {14, 12, 30, 12, 62}, {26, 24, 42, 24, 111}, {42, 40, 58, 60, 135}, {62, 60, 117, 90, 165}};

  std::cout << h.determinant() << '\n';

  Tensor::Matrix<int64_t> j = {{3, 4, 22, 4}, {14, 12, 30, 12}, {26, 24, 42, 24}, {42, 40, 58, 60}};

  std::cout << j.determinant() << '\n';

#endif

  Tensor::Matrix<double> k = {{-0.700624, 0.340395, -0.846329, -0.663135, -0.503968, 0.823883, 0.978628, -0.0576061},
                              {0.743513, 0.795462, -0.199831, 0.146875, 0.142725, -0.366635, 0.19976, 0.616893},
                              {0.713103, 0.489651, 0.836863, 0.755844, -0.421905, -0.785994, 0.894424, -0.0806125},
                              {-0.390112, -0.767406, 0.0110868, 0.244983, 0.24305, 0.309218, -0.36558, -0.794871},
                              {-0.0184016, 0.263631, 0.390246, 0.312225, -0.578232, -0.892475, 0.0635684, -0.346404},
                              {-0.383056, 0.11086, 0.357179, -0.872111, 0.794477, -0.334883, 0.701342, -0.953884},
                              {0.267538, 0.600714, -0.0252267, 0.589941, -0.223907, 0.786855, 0.329574, -0.863015},
                              {0.228132, -0.218749, 0.945731, -0.262868, 0.985404, 0.425773, -0.585217, 0.294802}};

  std::cout << "K DET: " << k.determinant() << '\n';

  Tensor::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};
  Tensor::Vector<std::complex<double>> cw = {2.0 - 3.0i, 4.0 + 1.0i, 2.0 + 5i};
  std::cout << cv << '\n';
  std::cout << cw << '\n';
  std::cout << Tensor::angle(cv, cw) << '\n';

  auto cvit = cv.begin();
  std::cout << *cvit << "\n";

  *cvit = 4.0 - 5.5i;
  for (auto& it : cv) {
    std::cout << it << "\n";
  }
#if 0
  std::cout << cv.norm() << '\n';
  std::cout << cv.rms() << '\n';
  std::cout << cv.mean() << '\n';
  std::cout << cv.stddev() << '\n';
  std::cout << cv.variance() << '\n';
  Tensor::Vector<int> v = {1, 2, 3, 4, 5};

  std::cout << v.stddev() << '\n';

  Tensor::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  std::cout << s * 4.3 << '\n';

  std::cout << s + s << '\n';

  Tensor::Matrix<long double> t = {{5, 7, 9}, {4, 3, 8}, {7, 5, 6}};
  auto ti = t.inverse();

  std::cout << "T INV:\n" << ti << '\n';

  std::cout << k.inverse() << "\n";
  // std::cout << t * ti << "\n";
  /*
    Tensor::Matrix<double> ss = {{1, 4}, {3, 8}};

    std::cout << ss.inverse() << "\n";*/

  Tensor::Vector<double> vf(4, gen);

  std::cout << vf << "\n";

  auto gl = [](const std::size_t i) -> double { return std::sqrt(3.0 * i + 7); };

  Tensor::Vector<double> vl(4, gl);

  std::cout << vl << "\n";

  auto mgl = [](const std::size_t i, const std::size_t j) -> double {
    return static_cast<double>(i * i + 1) / static_cast<double>(j + 1);
  };

  Tensor::Matrix<double> ml(3, mgl);

  std::cout << ml << "\n";

  DoubleFunctor f;

  Tensor::Matrix<double> mff(3, f);

  std::cout << mff << "\n";

  std::cout << mff(1, 2) << "\n";

  mff(2, 1) = 13;

  std::cout << mff << "\n";

  Tensor::Matrix<OrderedPair<int>> mal(
      3, [](const std::size_t i, const std::size_t j) { return OrderedPair<int>(std::make_pair(i, j)); });

  std::cout << "[" << mal[0][1].first << ", " << mal[0][1].second << "]\n";
  std::cout << mal[1][2] << "\n";

  std::cout << mal << "\n";
#endif
  // std::cout << mal << "\n";
  /* std::random_device r;
  std::default_random_engine eng(r());
  std::uniform_real_distribution<double> ud(0.1, 1.1);

  std::cout << ud(eng) << "\n";

  auto e = std::bind(ud, eng);

  std::cout << e() << "\n";

  Tensor::Matrix<double> mrd(500, std::bind(ud, eng));

  // std::cout << mrd << "\n";

  std::cout << mrd.determinant() << "\n";

  std::cout << mrd * mrd.inverse() << "\n"; */

  return 0;
}

double gen(const std::size_t i) { return i * i + 1; }