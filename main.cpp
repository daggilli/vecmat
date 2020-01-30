#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include "vecmat.h"

using namespace std::complex_literals;

int main() {
  /*
  VecMat::Vector<float> fv(3);

  VecMat::Matrix<double> dm(3);
  dm[2][0] = 1.8;

  std::cout << dm[0][1] << '\n';

  VecMat::Matrix<double> dom(3, 2);

  VecMat::Vector<double> dva = {3.0, 4, 7, 1.9}, dvb = {6.7, 1.1, 5.3, 6.1};
  std::cout << dva.size() << '\n';
  VecMat::Vector<double> dvc = dva + dvb;

  for (auto e : dvc) {
    std::cout << e << " ";
  }
  std::cout << '\n';

  VecMat::Matrix mil = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  std::cout << dvc << '\n';

  std::cout << mil << '\n';

  dvc += dvc;
  std::cout << dvc << '\n';

  VecMat::Matrix mila = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  VecMat::Matrix milb = {{10, 11, 12}, {13, 14, 15}, {16, 17, 18}};
  VecMat::Matrix milc = mila * milb;
  std::cout << milc << '\n';

  VecMat::Matrix a = {{1, 2, 3}, {4, 5, 6}};
  VecMat::Matrix b = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
  VecMat::Matrix c = a * b;
  std::cout << a << "\n\n";
  std::cout << b << "\n\n";
  std::cout << c << "\n"; */

  /*
    VecMat::Vector<double> v = {1, 2, 3, 4};
    VecMat::Vector<double> w = {1, 2, 3};
    VecMat::Matrix<double> n = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
    VecMat::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    VecMat::Vector<double> r = n * v;

    std::cout << r << "\n";

    r = w * n;

    std::cout << r << "\n";

    std::cout << w * w << '\n';

    std::cout << s * n << '\n';

    std::cout << VecMat::transpose(n) << '\n';

    auto x = s * n;

    std::cout << x << '\n';

    auto y = std::move(x);

    std::cout << y << '\n';
    std::cout << "mv\n";
    std::cout << x << '\n';

    VecMat::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};

    std::cout << cv << '\n';

    std::cout << cv.total() << '\n';

    std::cout << w.norm() << '\n';

    std::cout << cv.norm() << '\n';

    std::cout << w * 3.3 << '\n';
    std::cout << 6.2 * w << '\n';

    v = {1, 2, 3};
    w = {4, 5, 6};

    std::cout << VecMat::angle(v, w) << '\n';

    cv = {1.0 + 1.0i, 2.0 - 1.0i};
    VecMat::Vector<std::complex<double>> cw = {2.0 - 3.0i, 4.0 + 1.0i};

    std::cout << cv.norm() << " " << cw.norm() << '\n';

    std::cout << VecMat::angle(cv, cw) << '\n';

    VecMat::Matrix<std::complex<double>> cm = {{1.0 + 1.0i, 2.0 - 1.0i}, {2.0 - 3.0i, 4.0 + 1.0i}};

    std::cout << cm.trace() << '\n';

    */

#if 0
  VecMat::Vector<double> v = {1, 2, 3, 4};
  VecMat::Vector<double> w = {1, 2, 3};
  VecMat::Matrix<double> n = {{7, 8, 9, 10}, {11, 12, 13, 14}, {15, 16, 17, 18}};
  VecMat::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  // VecMat::Vector<double> r = n * v;
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

  VecMat::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};

  std::cout << cv << '\n';

  std::cout << cv.total() << '\n';

  w = v;

  std::cout << w << '\n';

  VecMat::Vector<double> vv(std::move(w));

  std::cout << vv << '\n';

  w = std::move(vv);

  std::cout << w << '\n';

  std::cout << w.norm() << '\n';

  cv = {1.0 + 1.0i, 2.0 - 1.0i};
  VecMat::Vector<std::complex<double>> cw = {2.0 - 3.0i, 4.0 + 1.0i};

  std::cout << cv.norm() << " " << cw.norm() << '\n';

  std::cout << VecMat::angle(cv, cw) << '\n';

  std::cout << std::boolalpha << s.isSquare() << '\n';
  std::cout << std::boolalpha << sxn.isSquare() << '\n';

  s[1][1] = 17;

  VecMat::Matrix<double> g = {{3, 4, 22}, {14, 12, 30}, {26, 24, 42}};
  std::cout << g.determinant() << '\n';

  VecMat::Matrix<double> f = {{3, 4, 22, 4, 54, 4, 153},          {14, 12, 30, 12, 62, 18, 165},
                              {26, 24, 42, 24, 111, 36, 183},     {42, 40, 58, 60, 135, 60, 207},
                              {62, 60, 117, 90, 165, 90, 237},    {86, 126, 153, 126, 201, 126, 273},
                              {171, 168, 195, 168, 243, 168, 315}};

  std::cout << f.determinant() << '\n';

  VecMat::Matrix<int64_t> h = {
      {3, 4, 22, 4, 54}, {14, 12, 30, 12, 62}, {26, 24, 42, 24, 111}, {42, 40, 58, 60, 135}, {62, 60, 117, 90, 165}};

  std::cout << h.determinant() << '\n';

  VecMat::Matrix<int64_t> j = {{3, 4, 22, 4}, {14, 12, 30, 12}, {26, 24, 42, 24}, {42, 40, 58, 60}};

  std::cout << j.determinant() << '\n';

#endif

  VecMat::Matrix<double> k = {{-0.700624, 0.340395, -0.846329, -0.663135, -0.503968, 0.823883, 0.978628, -0.0576061},
                              {0.743513, 0.795462, -0.199831, 0.146875, 0.142725, -0.366635, 0.19976, 0.616893},
                              {0.713103, 0.489651, 0.836863, 0.755844, -0.421905, -0.785994, 0.894424, -0.0806125},
                              {-0.390112, -0.767406, 0.0110868, 0.244983, 0.24305, 0.309218, -0.36558, -0.794871},
                              {-0.0184016, 0.263631, 0.390246, 0.312225, -0.578232, -0.892475, 0.0635684, -0.346404},
                              {-0.383056, 0.11086, 0.357179, -0.872111, 0.794477, -0.334883, 0.701342, -0.953884},
                              {0.267538, 0.600714, -0.0252267, 0.589941, -0.223907, 0.786855, 0.329574, -0.863015},
                              {0.228132, -0.218749, 0.945731, -0.262868, 0.985404, 0.425773, -0.585217, 0.294802}};

  std::cout << k.determinant() << '\n';

  VecMat::Vector<std::complex<double>> cv = {1.0 + 0.5i, 1.0i, 3.0 - 2.0i};

  std::cout << cv.norm() << '\n';
  std::cout << cv.rms() << '\n';
  std::cout << cv.mean() << '\n';
  std::cout << cv.stddev() << '\n';
  std::cout << cv.variance() << '\n';
  VecMat::Vector<int> v = {1, 2, 3, 4, 5};

  std::cout << v.stddev() << '\n';

  VecMat::Matrix<double> s = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  std::cout << s * 4.3 << '\n';

  std::cout << s + s << '\n';
  return 0;
}
