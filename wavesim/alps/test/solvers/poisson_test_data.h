#include <common/container/view_types.h>

template<typename T>
inline void get_solution(alps::MDView<T***>& d,
                         alps::MDView<T***>& du,
                         alps::MDView<T***>& dl,
                         alps::MDView<T***>& b,
                         alps::MDView<T***>& x)
{
  using namespace alps;

  std::vector<T> vec_d{1.96737369,
                       2.09413249,
                       1.98598118,
                       2.09547812,
                       1.90958619,
                       2.0530856,
                       1.99521273,
                       2.09563996,
                       1.90064632,
                       1.94602313,
                       1.95522112};
  std::vector<T> vec_du{-1.06779591,
                        -0.99794877,
                        -0.97270468,
                        -0.91582218,
                        -0.95632772,
                        -1.0609893,
                        -0.91330441,
                        -1.04160085,
                        -0.91496764,
                        -0.97324237,
                        0};
  std::vector<T> vec_dl{0,
                        -0.96746235,
                        -1.08435338,
                        -1.02573461,
                        -1.03540782,
                        -1.0836401,
                        -1.06939333,
                        -0.98648056,
                        -1.04860329,
                        -0.98793547,
                        -0.9045533};
  std::vector<T> vec_b{0.54701541,
                       0.10520587,
                       0.77325627,
                       0.02583852,
                       0.65943223,
                       0.89126474,
                       0.16106102,
                       0.66614267,
                       0.68407245,
                       0.69051037,
                       0.90028262};
  std::vector<T> vec_x{2.63758483,
                       4.34736596,
                       6.46024189,
                       7.54862368,
                       10.00810393,
                       11.12171139,
                       10.45946268,
                       9.65100791,
                       8.87176554,
                       6.62093081,
                       3.52352343};

  d            = MDView<T***>("", 4, 8, vec_d.size());
  du           = MDView<T***>("", 4, 8, vec_d.size());
  dl           = MDView<T***>("", 4, 8, vec_d.size());
  b            = MDView<T***>("", 4, 8, vec_d.size());
  x            = MDView<T***>("", 4, 8, vec_d.size());
  auto d_host  = create_mirror_view(d);
  auto du_host = create_mirror_view(du);
  auto dl_host = create_mirror_view(dl);
  auto b_host  = create_mirror_view(b);
  auto x_host  = create_mirror_view(x);

  for (int k = 0; k < d_host.extent_int(2); ++k) {
    for (int j = 0; j < d_host.extent_int(1); ++j) {
      for (int i = 0; i < d_host.extent_int(0); ++i) {
        d_host(i, j, k)  = (i + 1) * vec_d[k] / (j + 1);
        du_host(i, j, k) = (i + 1) * vec_du[k] / (j + 1);
        dl_host(i, j, k) = (i + 1) * vec_dl[k] / (j + 1);
        b_host(i, j, k)  = vec_b[k];
        x_host(i, j, k)  = (j + 1) * vec_x[k] / (i + 1);
      }
    }
  }

  deep_copy(d, d_host);
  deep_copy(du, du_host);
  deep_copy(dl, dl_host);
  deep_copy(x, x_host);
  deep_copy(b, b_host);
}

template<typename T>
inline void get_LU_result(alps::MDView<T***>& d,
                          alps::MDView<T***>& du,
                          alps::MDView<T***>& dl)
{
  using namespace alps;

  std::vector<T> vec_d{1.96737369,
                       1.56904041,
                       1.29630544,
                       1.32580083,
                       1.19435852,
                       1.18541057,
                       1.03806344,
                       1.22771897,
                       1.01100781,
                       1.05193609,
                       1.11833606};
  std::vector<T> vec_du{-1.06779591,
                        -0.99794877,
                        -0.97270468,
                        -0.91582218,
                        -0.95632772,
                        -1.0609893,
                        -0.91330441,
                        -1.04160085,
                        -0.91496764,
                        -0.97324237,
                        0};
  std::vector<T> vec_dl{0,
                        -0.49175322,
                        -0.69109334,
                        -0.7912754,
                        -0.78096785,
                        -0.90729884,
                        -0.90212907,
                        -0.95030855,
                        -0.85410694,
                        -0.97717888,
                        -0.85989378};

  d            = MDView<T***>("", 4, 8, vec_d.size());
  du           = MDView<T***>("", 4, 8, vec_d.size());
  dl           = MDView<T***>("", 4, 8, vec_d.size());
  auto d_host  = create_mirror_view(d);
  auto du_host = create_mirror_view(du);
  auto dl_host = create_mirror_view(dl);

  for (int k = 0; k < d_host.extent_int(2); ++k) {
    for (int j = 0; j < d_host.extent_int(1); ++j) {
      for (int i = 0; i < d_host.extent_int(0); ++i) {
        d_host(i, j, k)  = (i + 1) * vec_d[k] / (j + 1);
        dl_host(i, j, k) = vec_dl[k];
        du_host(i, j, k) = (i + 1) * vec_du[k] / (j + 1);
      }
    }
  }

  deep_copy(d, d_host);
  deep_copy(du, du_host);
  deep_copy(dl, dl_host);
}
