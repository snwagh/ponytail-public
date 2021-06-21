#include "seal_augment.h"

// #include "NTL/ZZ.h"
// #include "NTL/ZZ_p.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

// const int32_t LOG_N = 15;
// const int32_t LOG_P = 128;

using namespace std;

NTL::ZZ gen_p(int32_t log_p, int32_t log_n) {
  int32_t n = 1 << log_n;
  for (;;) {
    NTL::ZZ t = NTL::RandomBits_ZZ(log_p - log_n - 1);
    NTL::ZZ p = t * n * 2 + 1;
    if (NTL::ProbPrime(p) && NTL::NumBits(p) > log_p - 1) {
      return p;
    }
  }
}

int32_t rev(int32_t x, int32_t log_n) {
  int32_t ans = 0;
  for (int32_t i = 0; i < log_n; ++i) {
    ans |= ((x >> i) & 1) << (log_n - i - 1);
  }
  return ans;
}

// result[i] = sum_j a[j] omega^{ij}
// length of a must be 2^log_n

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

std::vector<NTL::ZZ_p> ntt_fast(const std::vector<NTL::ZZ_p> &a, int32_t log_n,
                                const NTL::ZZ_p &omega) {
  int32_t n = 1 << log_n;
  std::vector<NTL::ZZ_p> result(n);
  for (int32_t i = 0; i < n; ++i) {
    result[i] = a[rev(i, log_n)];
  }
  std::vector<NTL::ZZ_p> omega_two;
  NTL::ZZ_p cur(omega);
  for (int32_t i = 0; i < log_n; ++i) {
    omega_two.push_back(cur);
    cur *= cur;
  }
  for (int32_t i = 1; i <= log_n; ++i) {
    int32_t l = 1 << i;
    int32_t l2 = l / 2;
    NTL::ZZ_p ws = omega_two[log_n - i];
    NTL::ZZ_p wb = NTL::power(ws, l2);
    for (int32_t j = 0; j < n; j += l) {
      NTL::ZZ_p w1(1);
      NTL::ZZ_p w2(wb);
      for (int32_t k = j; k < j + l2; ++k) {
        NTL::ZZ_p tmp = result[k] + w1 * result[k + l2];
        result[k + l2] = result[k] + w2 * result[k + l2];
        result[k] = tmp;
        w1 *= ws;
        w2 *= ws;
      }
    }
  }
  return result;
}

void ntt_superfast(std::vector<NTL::ZZ_p> &a, int32_t log_n,
                   const std::vector<NTL::ZZ_p> &omega_pow, int32_t mul) {
  int32_t n = 1 << log_n;
  NTL::ZZ_p tmp;
  for (int32_t i = 0; i < n; ++i) {
    int32_t r = rev(i, log_n);
    if (r <= i) {
      continue;
    }
    tmp = a[i];
    a[i] = a[r];
    a[r] = tmp;
  }
  auto t1 = high_resolution_clock::now();
  for (int32_t i = 1; i <= log_n; ++i) {
    int32_t l = 1 << i;
    int32_t l2 = l / 2;
    for (int32_t j = 0; j < n; j += l) {
      for (int32_t k = j; k < j + l2; ++k) {
        const NTL::ZZ_p &w1 =
            omega_pow[(mul * (k - j) * (1 << (log_n - i)) + n) % n];
        const NTL::ZZ_p &w2 =
            omega_pow[(mul * (k - j + l2) * (1 << (log_n - i)) + n) % n];
        tmp = a[k] + w1 * a[k + l2];
        a[k + l2] = a[k] + w2 * a[k + l2];
        a[k] = tmp;
      }
    }
  }
  auto t2 = high_resolution_clock::now();
  std::cout << duration_cast<duration<double>>(t2 - t1).count() << std::endl;
}

std::vector<NTL::ZZ_p> ntt_inv_fast(const std::vector<NTL::ZZ_p> &a,
                                    int32_t log_n, const NTL::ZZ_p &omega) {
  NTL::ZZ_p omega_inv = NTL::power(omega, -1);
  int32_t n = 1 << log_n;
  NTL::ZZ_p n_inv = NTL::power(NTL::ZZ_p(n), NTL::to_ZZ(-1));
  auto b = ntt_fast(a, log_n, omega_inv);
  for (auto &x : b) {
    x *= n_inv;
  }
  return b;
}

void ntt_inv_superfast(std::vector<NTL::ZZ_p> &a, int32_t log_n,
                       const std::vector<NTL::ZZ_p> &omega_pow,
                       const NTL::ZZ_p &n_inv) {
  ntt_superfast(a, log_n, omega_pow, -1);
  for (auto &x : a) {
    x *= n_inv;
  }
}

// struct SEALStuff {
//   std::shared_ptr<seal::SEALContext> context;
//   std::unique_ptr<seal::Evaluator> evaluator;

// };

// struct SEALStuffClient {
//   std::unique_ptr<seal::KeyGenerator> key_generator;
//   seal::PublicKey public_key;
//   seal::SecretKey secret_key;
//   std::unique_ptr<seal::Encryptor> encryptor;
//   std::unique_ptr<seal::Decryptor> decryptor;
//   seal::GaloisKeys gal_keys;
// };

// struct BigPContext{
//     std::unique_ptr<SEALStuff> seal_stuff;
//     std::unique_ptr<SEALStuffClient> seal_stuff_client;
//     int32_t n;
//     NTL::ZZ p;
//     std::vector<NTL::ZZ> qs;
//     NTL::ZZ q;
//     NTL::ZZ_p xi;
//     NTL::ZZ_p xi_inv;
//     NTL::ZZ_p omega;
//     vector<int> index_map;

//     std::vector<uint64_t> mods1, mods2;
// };

// std::unique_ptr<SEALStuff> seal_stuff;
// std::unique_ptr<SEALStuffClient> seal_stuff_client;

std::unique_ptr<SEALStuff> setup_seal(size_t n, size_t logp,
                                      const std::vector<NTL::ZZ> &qs) {
  seal::EncryptionParameters encryption_parameters(seal::scheme_type::BFV);
  encryption_parameters.set_poly_modulus_degree(n);

  if (qs.size() > 0) {
    assert(logp > 60 && "Some issue with setup seal");
    std::vector<seal::SmallModulus> moduli;
    for (auto &x : qs) {
      std::ostringstream oss;
      oss << x;
      std::stringstream ss(oss.str());
      uint64_t val;
      ss >> val;
      std::cout << val << std::endl;
      moduli.push_back(seal::SmallModulus(val));
    }
    encryption_parameters.set_coeff_modulus(moduli);
  } else {
    encryption_parameters.set_coeff_modulus(seal::CoeffModulus::BFVDefault(n));
  }

  if (logp <= 60) {
    encryption_parameters.set_plain_modulus(
        seal::PlainModulus::Batching(n, logp));

    //  encryption_parameters.set_plain_modulus(1 << logp);
  } else {
    encryption_parameters.set_plain_modulus(
        seal::PlainModulus::Batching(n, 20));  // dummy
  }
  auto result = std::make_unique<SEALStuff>();
  result->context = seal::SEALContext::Create(encryption_parameters);
  result->parms =
      std::make_shared<seal::EncryptionParameters>(encryption_parameters);
  result->evaluator = std::make_unique<seal::Evaluator>(result->context);
  return result;
}

std::unique_ptr<SEALStuffClient> setup_seal_client(SEALStuff *seal_stuff,
                                                   int selection) {
  auto result = std::make_unique<SEALStuffClient>();
  result->key_generator =
      std::make_unique<seal::KeyGenerator>(seal_stuff->context);
  result->public_key = result->key_generator->public_key();
  result->secret_key = result->key_generator->secret_key();
  result->encryptor = std::make_unique<seal::Encryptor>(seal_stuff->context,
                                                        result->public_key);
  result->decryptor = std::make_unique<seal::Decryptor>(seal_stuff->context,
                                                        result->secret_key);

  seal::GaloisKeys gal_keys;

  /*****************************************************************
   Generate all the required Galois keys if GALKEY = 1;
   otherwise, generate power-of-two Galois keys as SEAL-default
   ******************************************************************/
  int n = seal_stuff->parms->poly_modulus_degree();
  int sub_optimal = n / (2 * YONGSOO_LIMIT * YONGSOO_LIMIT);

#if (GALKEY)
  if (selection == 1) {
    vector<int> steps(3 * (YONGSOO_LIMIT - 1), 0);
    for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
      steps[i] = (i + 1) * sub_optimal * YONGSOO_LIMIT;
      steps[i + (YONGSOO_LIMIT - 1)] = sub_optimal * (-YONGSOO_LIMIT + i + 1);
      steps[i + (YONGSOO_LIMIT - 1) * 2] = sub_optimal * (i + 1);
    }
    gal_keys = result->key_generator->galois_keys(steps);
  } else if ((selection == 2) || (selection == 3)) {
    vector<int> steps(YONGSOO_LIMIT - 1, 0);
    for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
      steps[i] = (i + 1) * sub_optimal * YONGSOO_LIMIT;
    }
    gal_keys = result->key_generator->galois_keys(steps);
  }
#else
#if (HOISTING)
  throw invalid_argument("hoisting should prepare all the rotation keys");
#endif
  gal_keys = result->key_generator->galois_keys();
#endif
  result->gal_keys = gal_keys;
  // result->gal_keys = result->key_generator->galois_keys();

  result->relin_keys = result->key_generator->relin_keys();

  return result;
}

uint64_t to_uint64_t(const NTL::ZZ &x) {
  return NTL::conv<uint64_t>(x);
  // std::ostringstream oss;
  // oss << x;
  // std::stringstream ss(oss.str());
  // uint64_t res;
  // ss >> res;
  // return res;
}

NTL::ZZ from_uint64_t(const uint64_t &x) {
  std::ostringstream oss;
  oss << x;
  NTL::ZZ res = NTL::to_ZZ(oss.str().c_str());
  return res;
}

NTL::ZZ_p to_zzp(const NTL::ZZ &x) {
  std::ostringstream oss;
  oss << x;
  std::stringstream ss(oss.str());
  NTL::ZZ_p res;
  ss >> res;
  return res;
}

void test_ntt(int logp, int logn) {
  cout << "Benchmarking NTT...." << endl;
  int32_t n = 1 << logn;
  auto p = gen_p(logp, logn);
  std::cout << p << std::endl;

  NTL::ZZ_p::init(p);
  NTL::ZZ_p xi;
  for (int32_t i = 2;; ++i) {
    xi = NTL::power(NTL::ZZ_p(i), (p - 1) / (2 * n));
    if (NTL::power(xi, n) != 1) {
      break;
    }
  }

  if (NTL::power(xi, 2 * n) != 1) {
    cout << "WARNING: somethign wrong" << endl;
  }

  std::cout << "xi = " << xi << std::endl;
  NTL::ZZ_p omega = xi * xi;

  NTL::ZZ_p xi_inv = NTL::power(xi, -1);
  std::vector<NTL::ZZ_p> a;  // original vector

  NTL::ZZ_p cur(1);
  for (int32_t i = 0; i < n; ++i) {
    a.push_back(NTL::random_ZZ_p());
  }
  std::vector<NTL::ZZ_p> a_tilde(n);
  for (int32_t i = 0; i < n; i++) {
    a_tilde[i] = a[i] * cur;
    cur *= xi;
    // cout << a[i] << ", ";
  }
  std::cout << std::endl;
  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);
  for (int i = 0; i < 100; i++) auto b = ntt_fast(a_tilde, logn, omega);
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "100 ntt time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
       << endl;

  // inverse
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total_inv(0);
  for (int i = 0; i < 100; i++) auto b = ntt_inv_fast(a_tilde, logn, omega);
  time_end = chrono::high_resolution_clock::now();
  time_total_inv =
      chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "100 inverse ntt time: \t\t" << time_total_inv.count()
       << PRINT_STRING(TIME_UNIT) << endl;

  // std::vector<NTL::ZZ_p> bcheck(n);
  // for (int32_t i = 0; i < n; ++i) {
  //   bcheck[i] = 0;
  //   for (int32_t j = 0; j < n; j++) {
  //     bcheck[i] += a[j] * NTL::power(xi, (2 * i + 1) * j);
  //   }
  //   std::cout << "(should be zero: " << i << ", " << b[i] - bcheck[i] << ")
  //   ";
  // }
  // std::cout << std::endl;
}

seal::Ciphertext big_encrypt(const std::vector<NTL::ZZ_p> &plaintext,
                             std::shared_ptr<BigPContext> context) {
  int n = context->n;

  std::vector<NTL::ZZ_p> plain_permuted(n);  // original vector

  for (int32_t i = 0; i < n; ++i) {
    plain_permuted[context->index_map[i]] =
        plaintext[i];  //.push_back(NTL::random_ZZ_p());
  }
  int logn = log2(n);

  auto b =
      ntt_inv_fast(plain_permuted, logn, context->omega);  // inverse ntt of a
  NTL::ZZ_p cur(1);
  for (int32_t i = 0; i < n; ++i) {
    b[i] *= cur;
    cur *= context->xi_inv;
  }
  // auto b = plaintext;
  // seal::Plaintext zero(n);
  seal::Ciphertext zero_ciphertext;
  // context->seal_stuff_client->encryptor->encrypt(zero, zero_ciphertext);
  context->seal_stuff_client->encryptor->encrypt_zero(zero_ciphertext);

  vector<NTL::ZZ> scaled_b(n);
  for (int j = 0; j < n; j++) {
    NTL::ZZ bj = NTL::rep(b[j]);
    scaled_b[j] =
        (bj * context->q) /
        context
            ->p;  // should be rounding instead of flooring. but we don't care.
  }

  for (size_t i = 0; i < context->qs.size() - 1; ++i) {
    // auto delta = (q /  p) % qs[i];
    auto qsill = NTL::conv<uint64_t>(context->qs[i]);
    for (int32_t j = 0; j < n; ++j) {
      auto cur = scaled_b[j] % context->qs[i];
      zero_ciphertext.data_[i * n + j] += NTL::conv<uint64_t>(cur);
      zero_ciphertext.data_[i * n + j] %= qsill;
    }
  }
  return zero_ciphertext;
}

seal::Plaintext big_encode(const std::vector<NTL::ZZ_p> &plaintext,
                           std::shared_ptr<BigPContext> context) {
  int n = context->n;
  NTL::ZZ_p omega = context->omega;
  NTL::ZZ_p xi_inv = context->xi_inv;

  // todo: implement this.
  std::vector<NTL::ZZ_p> plain_permuted(n);  // original vector

  for (int32_t i = 0; i < n; ++i) {
    plain_permuted[context->index_map[i]] =
        plaintext[i];  //.push_back(NTL::random_ZZ_p());
  }
  int logn = log2(n);
  auto b = ntt_inv_fast(plain_permuted, logn, omega);  // inverse ntt of a
  NTL::ZZ_p cur(1);
  for (int32_t i = 0; i < n; ++i) {
    b[i] *= cur;
    cur *= xi_inv;
  }
  // auto b = plaintext;
  seal::Plaintext result(0);
  size_t coeff_mod_count = context->qs.size() - 1;
  result.resize(n * coeff_mod_count);
  // cout << "size of plaintext = " << n << ", " << coeff_mod_count << endl;
  for (size_t i = 0; i < coeff_mod_count; i++) {
    auto qi = context->seal_stuff->parms->coeff_modulus()[i].value();
    for (int j = 0; j < n; j++) {
      result.data()[i * n + j] =
          NTL::conv<uint64_t>(NTL::conv<NTL::ZZ>(b[j]) % context->qs[i]);
    }
  }
  return result;
}

std::vector<NTL::ZZ_p> big_decode(const seal::Plaintext &plaintext,
                                  std::shared_ptr<BigPContext> context) {
  // fixme: implement this.
  int n = context->n;
  NTL::ZZ p = context->p;
  NTL::ZZ q = context->q;
  NTL::ZZ noise(0);

  vector<NTL::ZZ_p> temp(n), result(n);

  for (int32_t i = 0; i < n; ++i) {
    auto cur_value = NTL::conv<NTL::ZZ>(*plaintext.data(i));
    auto cur_modulus = context->qs[0];
    for (size_t j = 1; j < context->qs.size() - 1; ++j) {
      NTL::CRT(cur_value, cur_modulus,
               NTL::conv<NTL::ZZ>(*plaintext.data(j * n + i)), context->qs[j]);
    }
    cur_value += cur_modulus;
    if (cur_value < 0) {
      std::cerr << "Invalid range" << std::endl;
      exit(1);
    }
    while (cur_value >= cur_modulus) {
      cur_value -= cur_modulus;
    }
    temp[i] = NTL::conv<NTL::ZZ_p>(cur_value);
  }

  int logn = log2(n);

  NTL::ZZ_p cur1(1);
  for (int32_t i = 0; i < n; ++i) {
    temp[i] *= cur1;
    cur1 *= context->xi;
  }

  auto temp2 = ntt_fast(temp, logn, context->omega);  // decryption
  // std::vector<NTL::ZZ_p> aa_permuted(n);  // original vector

  for (int32_t i = 0; i < n; ++i) {
    result[i] = temp2[context->index_map[i]];
  }
  return result;
}

std::vector<NTL::ZZ_p> big_decrypt(const seal::Ciphertext &ciphertext,
                                   std::shared_ptr<BigPContext> context) {
  std::vector<std::vector<uint64_t>> result;
  context->seal_stuff_client->decryptor->decrypt_hacked(ciphertext, result);
  std::vector<NTL::ZZ_p> bb;

  int n = context->n;
  NTL::ZZ p = context->p;
  NTL::ZZ q = context->q;
  NTL::ZZ noise(0);
  for (int32_t i = 0; i < n; ++i) {
    auto cur_value = NTL::conv<NTL::ZZ>(result[0][i]);
    auto cur_modulus = context->qs[0];
    for (size_t j = 1; j < context->qs.size() - 1; ++j) {
      NTL::CRT(cur_value, cur_modulus, NTL::conv<NTL::ZZ>(result[j][i]),
               context->qs[j]);
    }
    cur_value += cur_modulus;
    if (cur_value < 0) {
      std::cerr << "Invalid range" << std::endl;
      exit(1);
    }
    while (cur_value >= cur_modulus) {
      cur_value -= cur_modulus;
    }
    // cur_value is in [0, q-1].
    // compute round(t*cur_value / q)

    auto di = p * cur_value / q;
    auto mo = (p * cur_value) % q;
    if (mo > q / 2) {
      di += 1;
      mo -= q;
    }

    // auto di = cur_value / delta;
    // auto mo = cur_value % delta;
    // if (mo > delta2) {
    //   di += 1;
    //   mo -= delta;
    // }
    noise = max(noise, mo);
    bb.push_back(to_zzp(di));
  }
  cout
      << "noise budget (needs to be positive integer for correct decryption =  "
      << NTL::NumBits(q) - NTL::NumBits(noise) - 1 << endl;
  NTL::ZZ_p cur1(1);
  for (int32_t i = 0; i < n; ++i) {
    bb[i] *= cur1;
    cur1 *= context->xi;
  }
  int logn = log2(n);
  auto aa = ntt_fast(bb, logn, context->omega);  // decryption
  // std::vector<NTL::ZZ_p> aa_permuted(n);  // original vector

  for (int32_t i = 0; i < n; ++i) {
    bb[i] = aa[context->index_map[i]];
  }

  return bb;
}

// compute noise budget.

bool check_equal_zzp(vector<NTL::ZZ_p> &a, vector<NTL::ZZ_p> &b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

void pt_ct_mul(seal::Ciphertext &ct, const seal::Plaintext &pt,
               std::shared_ptr<BigPContext> context) {
  context->seal_stuff->evaluator->multiply_plain_inplace_hacked(ct, pt);
}

void test_pt_ct(std::shared_ptr<BigPContext> context) {
  int n = context->n;

  std::vector<NTL::ZZ_p> a;
  for (int32_t i = 0; i < n; ++i) {
    a.push_back(NTL::random_ZZ_p());
  }

  std::vector<NTL::ZZ_p> b;
  for (int32_t i = 0; i < n; ++i) {
    b.push_back(NTL::random_ZZ_p());
  }

  auto a_ct = big_encrypt(a, context);
  auto dummy = big_decrypt(a_ct, context);

  auto b_pt = big_encode(b, context);  // without ntt
  // encrypt and then decrypt

  context->seal_stuff->evaluator->multiply_plain_inplace_hacked(a_ct, b_pt);

  auto result = big_decrypt(a_ct, context);

  std::vector<NTL::ZZ_p> check(n);
  for (int i = 0; i < n; i++) {
    check[i] = a[i] * b[i];
  }

  cout << "ok ? " << check_equal_zzp(result, check) << endl;
}

void test_encrypt_decrypt() {
  cout << "Testing encrypt / decrypt ... ";
  std::shared_ptr<BigPContext> context = setup_params(128, 15, 12, 60, 2);
  // test encryption and decryption
  int n = context->n;
  int iterations = 10;

  std::vector<NTL::ZZ_p> a;
  for (int32_t i = 0; i < n; ++i) {
    a.push_back(NTL::random_ZZ_p());
  }

  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);
  //for (int i = 0; i < 100; i++) plain = big_encode(test_vec, bigpcontext);
  seal::Ciphertext a_ct;
  for (size_t i = 0 ; i < iterations; i++) 
    a_ct = big_encrypt(a, context); 
  time_end = chrono::high_resolution_clock::now();
    time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << iterations << " encrypt  time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;
  
  // encrypt and then decrypt
  time_start = chrono::high_resolution_clock::now();
  vector<NTL::ZZ_p> a_dec;
  for (size_t i = 0 ; i < iterations; i++) 
    a_dec = big_decrypt(a_ct, context);
  time_end = chrono::high_resolution_clock::now();
    time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << iterations << " decrypt time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;
  // cout << "ok ? " << check_equal_zzp(a, a_dec) << endl;
}

std::shared_ptr<BigPContext> setup_params(size_t logp, size_t logn,
                                          size_t num_q, size_t logqi,
                                          int selection) {
  auto result = std::make_shared<BigPContext>();
  int n = 1 << logn;
  result->n = n;
  result->p = gen_p(logp, logn);
  std::cout << "log p = " << NTL::NumBits(result->p) << std::endl;
  result->qs.clear();
  for (int32_t i = 0; i < num_q; ++i) {
    result->qs.push_back(gen_p(logqi, logn));
  }
  result->q = 1;
  for (int32_t i = 0; i < num_q - 1; ++i) {
    result->q *= result->qs[i];
  }
  std::cout << " log q = " << NTL::NumBits(result->q) << endl;

  NTL::ZZ_p::init(result->p);
  NTL::ZZ_p xi, omega, xi_inv;
  for (int32_t i = 2;; ++i) {
    xi = NTL::power(NTL::ZZ_p(i), (result->p - 1) / (2 * n));
    if (NTL::power(xi, n) != 1) {
      break;
    }
  }

  if (NTL::power(xi, 2 * n) != 1) {
    cout << "WARNING: somethign wrong" << endl;
  }

  omega = xi * xi;
  xi_inv = NTL::power(xi, -1);

  result->xi = xi;
  result->omega = omega;
  result->xi_inv = xi_inv;

  result->index_map.resize(n);
  int gen = 3;
  int pos = 1;
  int m = 2 * n;
  for (int i = 0; i < n / 2; i++) {
    result->index_map[i] = (pos - 1) >> 1;
    result->index_map[n / 2 + i] = (m - pos - 1) >> 1;
    pos *= gen;
    pos %= m;
  }
  result->seal_stuff = setup_seal(n, logp, result->qs);

  result->seal_stuff_client =
      setup_seal_client(result->seal_stuff.get(), selection);

  auto context_data =
      result->seal_stuff->evaluator->context_->first_context_data();
  auto &parms = context_data->parms();
  auto &coeff_modulus = parms.coeff_modulus();
  auto &base_converter = context_data->base_converter();
  auto &bsk_modulus = base_converter->get_bsk_mod_array();
  size_t coeff_mod_count = coeff_modulus.size();
  size_t bsk_base_mod_count = base_converter->bsk_base_mod_count();
  result->mods1.resize(coeff_mod_count);
  result->mods2.resize(bsk_base_mod_count);
  for (size_t i = 0; i < coeff_mod_count; ++i) {
    result->mods1[i] =
        to_uint64_t(result->p % from_uint64_t(coeff_modulus[i].value()));
    // std::cout << mods1[i] << " ";
  }
  // std::cout << std::endl;
  for (size_t i = 0; i < bsk_base_mod_count; ++i) {
    result->mods2[i] =
        to_uint64_t(result->p % from_uint64_t(bsk_modulus[i].value()));
    // std::cout << mods2[i] << " ";
  }
  print_parameters(result->seal_stuff->context);
  cout << "Actual plain modulus " << result->p << endl;

  result->omega_pow.resize(n);
  result->omega_pow[0] = 1;
  for (int32_t i = 1; i < n; ++i) {
    result->omega_pow[i] = result->omega_pow[i - 1] * result->omega;
  }
  result->n_inv = NTL::power(NTL::ZZ_p(n), NTL::to_ZZ(-1));

  return result;
}

void rotate(seal::Ciphertext &ciphertext, int32_t step,
            std::shared_ptr<BigPContext> context) {
  context->seal_stuff->evaluator->rotate_rows_inplace(
      ciphertext, step, context->seal_stuff_client->gal_keys);
}

void test_rotate(std::shared_ptr<BigPContext> context) {
  NTL::ZZ_p::init(context->p);
  int n = context->n;
  std::vector<NTL::ZZ_p> a;
  for (int32_t i = 0; i < n; ++i) {
    a.push_back(NTL::random_ZZ_p());
  }
  auto a_ct = big_encrypt(a, context);
  // auto b_ct = big_encrypt(b);

  int step = 3;

  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);
  //for (int i = 0; i < 100; i++) plain = big_encode(test_vec, bigpcontext);
  
  for (size_t i = 0 ; i < 100; i++) 
    context->seal_stuff->evaluator->rotate_rows_inplace(
        a_ct, step, context->seal_stuff_client->gal_keys);
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "100 rotate (without hoisting) time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;

  auto adec = big_decrypt(a_ct, context);

  for (int32_t i = 0; i < n / 2; ++i) {
    if (adec[i] != a[(i + step) % (n / 2)]) {
      std::cerr << "Invalid rotation" << std::endl;
      exit(1);
    }
  }
  return;
}

void ct_ct_mul(const seal::Ciphertext &a, const seal::Ciphertext &b,
               seal::Ciphertext &c, std::shared_ptr<BigPContext> context) {
  context->seal_stuff->evaluator->multiply_hacked(a, b, c, context->mods1,
                                                  context->mods2);
}

void test_multiply() {
  cout << "Testing ct multiply...";
  std::shared_ptr<BigPContext> context = setup_params(128, 15, 12, 60, 2);
  int n = context->n;
  NTL::ZZ p = context->p;
  NTL::SetSeed((NTL::conv<NTL::ZZ>((long)time(NULL))));
  std::vector<NTL::ZZ_p> a;
  for (int32_t i = 0; i < n; ++i) {
    a.push_back(NTL::random_ZZ_p());
  }
  // cout << "a0 = " << a[0] << endl;
  std::vector<NTL::ZZ_p> b;
  for (int32_t i = 0; i < n; ++i) {
    b.push_back(NTL::random_ZZ_p());
  }
  auto a_ct = big_encrypt(a, context);
  auto b_ct = big_encrypt(b, context);

  // just for checking noise budget in inputs.
  auto aa = big_decrypt(a_ct, context);
  auto bb = big_decrypt(b_ct, context);

  // auto context_data =
  // context->seal_stuff->evaluator->context_->get_context_data(a_ct.parms_id());
  // auto &parms = context_data->parms();
  // auto &coeff_modulus = parms.coeff_modulus();
  // auto &base_converter = context_data->base_converter();
  // auto &bsk_modulus = base_converter->get_bsk_mod_array();
  // size_t coeff_mod_count = coeff_modulus.size();
  // size_t bsk_base_mod_count = base_converter->bsk_base_mod_count();
  // std::vector<uint64_t> mods1(coeff_mod_count);
  // std::vector<uint64_t> mods2(bsk_base_mod_count);
  // for (size_t i = 0; i < coeff_mod_count; ++i) {
  //   mods1[i] = to_uint64_t(p % from_uint64_t(coeff_modulus[i].value()));
  //   //std::cout << mods1[i] << " ";
  // }
  // //std::cout << std::endl;
  // for (size_t i = 0; i < bsk_base_mod_count; ++i) {
  //   mods2[i] = to_uint64_t(p % from_uint64_t(bsk_modulus[i].value()));
  //   //std::cout << mods2[i] << " ";
  // }
  // std::cout << std::endl;
  seal::Ciphertext c_ct;
  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);
  for (int i = 0; i < 100; i++) 
    context->seal_stuff->evaluator->multiply_hacked(
      a_ct, b_ct, c_ct, context->mods1, context->mods2);
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "100 multiply time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;

  auto c = big_decrypt(c_ct, context);
  for (int32_t i = 0; i < n; ++i) {
    if (a[i] * b[i] - c[i] != 0) {
      cout << "multiply result wrong" << std::endl;
      break;
    }
  }
}

/*
int main() {

  std::shared_ptr<BigPContext> context = setup_params();

  //test_pt_ct(context);

  // Hao: first test if encryption and decryption work.
  //test_encrypt_decrypt(context);

  //test_rotate(context);



  for(int i = 0; i < 2; i++){
      cout << i << "-th trial" << endl;
      test_multiply(context);
  }

  return 0;
  n = 1 << LOG_N;
  p = gen_p(LOG_P, LOG_N);
  std::cout << p << std::endl;
  for (int32_t i = 0; i < NUM_Q; ++i) {
    qs.push_back(gen_p(LOG_QI, LOG_N));
  }
  q = 1;
  for (int32_t i = 0; i < NUM_Q - 1; ++i) {
    q *= qs[i];
  }
  std::cout << "log q = " << NTL::NumBits(q) << std::endl;
  NTL::ZZ_p::init(p);
  for (int32_t i = 2;; ++i) {
    xi = NTL::power(NTL::ZZ_p(i), (p - 1) / (2 * n));
    if (NTL::power(xi, n) != 1) {
      break;
    }
  }
  seal_stuff = setup_seal(qs);
  seal_stuff_client = setup_seal_client(seal_stuff.get());

  std::cout << "log p = " << NTL::NumBits(p) << std::endl;
  vector<int> index_map(n);
  int gen = 3;
  int pos = 1;
  int m = 2 * n;
  for (int i = 0; i < n / 2; i++) {
    index_map[i] = (pos - 1) >> 1;
    index_map[n / 2 + i] = (m - pos - 1) >> 1;
    pos *= gen;
    pos %= m;
  }

  omega = xi * xi;
  xi_inv = NTL::power(xi, -1);

  for (;;) {

    break;
  }

  return 0;
}
*/

std::shared_ptr<SmallPContext> setup_params(size_t logp, size_t logn,
                                            int selection) {
  auto result = std::make_shared<SmallPContext>();

  size_t n = 1 << logn;

  result->seal_stuff = setup_seal(n, logp, vector<NTL::ZZ>());

  result->seal_stuff_client =
      setup_seal_client(result->seal_stuff.get(), selection);

  // EncryptionParameters parms(scheme_type::BFV);
  // parms.set_poly_modulus_degree(poly_modulus_degree);
  // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  // parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, logp));
  // shared_ptr<SEALContext> context = SEALContext::Create(parms);

  // KeyGenerator keygen(context);
  // auto secret_key = keygen.secret_key();
  // auto public_key = keygen.public_key();

  // Encryptor encryptor(context, public_key);
  // Decryptor decryptor(context, secret_key);
  // Evaluator evaluator(context);
  //
  result->batch_encoder =
      std::make_unique<seal::BatchEncoder>(result->seal_stuff->context);
  // size_t slot_count = batch_encoder.slot_count();
  // RelinKeys relin_keys = keygen.relin_keys();
  print_parameters(result->seal_stuff->context);
  // cout << endl;

  return result;
}

void test_encode_decode() {
  cout << "Testing encode/decode for big p ..." << endl;
  std::shared_ptr<BigPContext> bigpcontext = setup_params(128, 15, 12, 60, 2);

  vector<NTL::ZZ_p> test_vec(
      bigpcontext->seal_stuff->parms->poly_modulus_degree());
  for (size_t i = 0; i < test_vec.size(); i++) {
    test_vec[i] = NTL::random_ZZ_p();
  }

  seal::Plaintext plain;
  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);
  for (int i = 0; i < 100; i++) plain = big_encode(test_vec, bigpcontext);
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "100 big_encode time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;

  // seal::Plaintext plain = big_encode(test_vec, bigpcontext);

  auto vec = big_decode(plain, bigpcontext);

  cout << "ok ? " << check_equal_zzp(vec, test_vec) << endl;
}



void test_linear_monomials(size_t u, size_t v){
  cout << "Testing linear monomials for big p ..." << endl;
  std::shared_ptr<BigPContext> bigpcontext = setup_params(128, 15, 12, 60, 2);


  int n = bigpcontext->seal_stuff->parms->poly_modulus_degree(); 

  vector<NTL::ZZ_p> test_vec(n);
  for (size_t i = 0; i < test_vec.size(); i++) {
    test_vec[i] = NTL::random_ZZ_p();
  }
  vector<seal::Ciphertext> cts_u(u); 
  vector<seal::Ciphertext> cts_v(v); 

  for (size_t i = 0; i < u; i++){
    cts_u[i] = big_encrypt(test_vec, bigpcontext); 
  }
  for (size_t j = 0; j < v; j++){
    cts_v[j] = big_encrypt(test_vec, bigpcontext); 
  }
  chrono::high_resolution_clock::time_point time_start, time_end;
  time_start = chrono::high_resolution_clock::now();
  chrono::TIME_UNIT time_total(0);

  random_device rd; 

  for (size_t j = 0; j < u; j ++){
    for(size_t i  = 0; i < v; i++) {
        bigpcontext->seal_stuff->evaluator->multiply_monomial(cts_u[j],  rd() % (2*n));
        bigpcontext->seal_stuff->evaluator->add_inplace(cts_v[i], cts_u[j]);
    }
  }

  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "u, v verification time \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;

  // seal::Plaintext plain = big_encode(test_vec, bigpcontext);

  // auto vec = big_decode(plain, bigpcontext);

  // cout << "ok ? " << check_equal_zzp(vec, test_vec) << endl;
}