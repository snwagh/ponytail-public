#include "seal_augment.h"

#include <chrono>
#include <cstdint>

const int32_t log_n = 15;
const int32_t n = 1 << log_n;

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

int32_t main() {
  auto params = setup_params(128, log_n, 12, 60, 1);
   
  //test_rotate(params); 
  //return 0;

  NTL::ZZ_p::init(params->p);
  std::vector<NTL::ZZ_p> a(n);
  for (int32_t i = 0; i < n; ++i) {
    a[i] = NTL::random_ZZ_p();
  }
  auto t1 = high_resolution_clock::now();
  auto b = ntt_fast(a, log_n, params->omega);
  auto t2 = high_resolution_clock::now();
  auto c = ntt_inv_fast(b, log_n, params->omega);
  auto t3 = high_resolution_clock::now();
  for (int32_t i = 0; i < n; ++i) {
    if (c[i] != a[i]) {
      std::cerr << "Wrong" << std::endl;
      exit(1);
    }
  }
  auto aa = a;
  auto t4 = high_resolution_clock::now();
  ntt_superfast(aa, log_n, params->omega_pow, 1);
  auto t5 = high_resolution_clock::now();
  for (int32_t i = 0; i < n; ++i) {
    if (aa[i] != b[i]) {
      std::cerr << "Wrong 2" << std::endl;
      exit(1);
    }
  }
  auto t6 = high_resolution_clock::now();
  ntt_inv_superfast(aa, log_n, params->omega_pow, params->n_inv);
  auto t7 = high_resolution_clock::now();
  for (int32_t i = 0; i < n; ++i) {
    if (aa[i] != c[i]) {
      std::cerr << "Wrong 3" << std::endl;
      exit(1);
    }
  }
  std::cout << "OK" << std::endl;
  std::cout << "Old NTT: " << duration_cast<duration<double>>(t2 - t1).count()
            << std::endl;
  std::cout << "Old inverse NTT: "
            << duration_cast<duration<double>>(t3 - t2).count() << std::endl;
  std::cout << "New NTT: " << duration_cast<duration<double>>(t5 - t4).count()
            << std::endl;
  std::cout << "New inverse NTT: "
            << duration_cast<duration<double>>(t7 - t6).count() << std::endl;
  return 0;
}
