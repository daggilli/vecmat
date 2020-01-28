#include <complex>
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
  return 0;
}
