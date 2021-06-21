/*******************************************************
 * 2019 Sameer Wagh <snwagh@gmail.com>
 *
 * Ciphertext-ciphertext and Plaintext-Ciphertext matrix
 * multiplication implemented using SEAL and BFV in this
 * file. This assumes columnwise encoding of matrices.
 * Algorithm used is Yongsoo's CCS paper.
 *
 * This file is part of ponytail project.
 *******************************************************/

#pragma once

#include <assert.h>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include "seal/seal.h"
#include "seal_augment.h"

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"

#define POLY_MODULUS_DEGREE (32768)
// #define TIME_UNIT milliseconds
//#define STRING(X) #X
// #define PRINT_STRING(X) " " STRING(X)
#define VERBOSE false

#define NUM_Q (12)
#define LOG_QI (60)

/* optimiation options:
 * 1. GALKEY: generate all the required key if true; otherwise, generate
 * power-of-two keys as SEAL-default.
 * 2. HOISTING: Apply hoisting technique to generate encryptions of permuted
 * matrices if true
 * 3. FASTBLOCK: Instead of computing one-by-one subblocks, we can perform an
 * expensive pre-computation that depends on input subslocks.
 */
#define GALKEY true
#define FASTBLOCK true
#define HOISTING true

/*
 * Single ciphertext Yongsoo matrix multiplication. Functon is overloaded with
 * two types, one with parameter A as vector of plaintexts and other with
 * parameter A as a vector of ciphertexts.
 *
 * @param A:             Plaintext/Ciphertext A representing a square matrix of
 *                       dimension YONGSOO_LIMIT. Matrix A is values of
 * \tilde{A} i.e., the output of the first transformation completed.
 * @param B:             Ciphertext B representing the second square matrix of
 * same dimension. Note that B has to be repeated ciphertext.
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 * @param keygen:        HE KeyGenerator
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix
 */
void matrix_mult_unit(const seal::Plaintext& A, const seal::Ciphertext& B,
                      seal::Ciphertext& C, seal::Evaluator& evaluator,
                      seal::BatchEncoder& batch_encoder,
                      seal::GaloisKeys& gal_keys);
void matrix_mult_unit(const seal::Ciphertext& A, const seal::Ciphertext& B,
                      seal::Ciphertext& C, seal::Evaluator& evaluator,
                      seal::BatchEncoder& batch_encoder,
                      seal::GaloisKeys& gal_keys);

void matrix_mult_unit_colwise(
    const std::vector<std::vector<seal::Ciphertext> >& A_cts,
    const std::vector<std::vector<seal::Ciphertext> >& B_col_cts,
    std::vector<seal::Ciphertext>& C, seal::Evaluator& evaluator,
    seal::RelinKeys& relin_keys,
    std::shared_ptr<BigPContext> bigpcontext = nullptr);

void matrix_mult_unit_colwise(
    const std::vector<std::vector<seal::Plaintext> >& A_pts,
    const std::vector<std::vector<seal::Ciphertext> >& B_col_cts,
    std::vector<seal::Ciphertext>& C, seal::Evaluator& evaluator);

void matrix_mult_unit_rowwise(
    const std::vector<std::vector<seal::Ciphertext> >& A_row_cts,
    const std::vector<std::vector<seal::Ciphertext> >& B_cts,
    std::vector<seal::Ciphertext>& C, seal::Evaluator& evaluator,
    seal::RelinKeys& relin_keys);

void matrix_mult_unit_rowwise(
    const std::vector<std::vector<seal::Plaintext> >& A_row_pts,
    const std::vector<std::vector<seal::Ciphertext> >& B_cts,
    std::vector<seal::Ciphertext>& C, seal::Evaluator& evaluator,
    std::shared_ptr<BigPContext> bigpcontext = nullptr);

/*
 * Extending Yongsoo matrix multiplication to large blocks. Overloaded with two
 * types, one with parameter A as vector of plaintexts and other with parameter
 * A as a vector of ciphertexts.
 *
 * @param A:               vector of plaintext A representing matrix of size
 * k1*k2. Each plaintext has YONGSOO_LIMIT number of data values. Matrix A is
 * values of \tilde{A} i.e., the output of the first transformation completed.
 * All such matrices are in row major blocks.
 * @param B:               vector of ciphertext B representing matrix of size
 * k2*k3
 * @param k1, k2, k3:      Matrix sizes
 * @param evaluator:       HE Evaluator
 * @param batch_encoder:   HE batch_encoder
 * @param encrypted_zeros: Chphertext containing zeros.
 *
 * @return C:              vector of ciphertext C=AxB representing the product
 * matrix size of C is k1*k3
 */
void matrix_mult_general(const std::vector<seal::Plaintext>& A,
                         const std::vector<seal::Ciphertext>& B,
                         std::vector<seal::Ciphertext>& C, const size_t k1,
                         const size_t k2, const size_t k3,
                         seal::Evaluator& evaluator,
                         seal::BatchEncoder& batch_encoder,
                         seal::GaloisKeys& gal_keys);

void matrix_mult_general_ntt(const std::vector<seal::Plaintext>& A,
                             const std::vector<seal::Ciphertext>& B,
                             std::vector<seal::Ciphertext>& C, const size_t k1,
                             const size_t k2, const size_t k3,
                             seal::Evaluator& evaluator,
                             seal::BatchEncoder& batch_encoder,
                             seal::GaloisKeys& gal_keys);

void matrix_mult_general(const std::vector<seal::Ciphertext>& A,
                         const std::vector<seal::Ciphertext>& B,
                         std::vector<seal::Ciphertext>& C, const size_t k1,
                         const size_t k2, const size_t k3,
                         seal::Evaluator& evaluator,
                         seal::BatchEncoder& batch_encoder,
                         seal::GaloisKeys& gal_keys,
                         seal::RelinKeys& relin_keys);

/*
 * Fast version of blockwise matrix computation
 */
void matrix_mult_general_fast_ptct(
    const std::vector<seal::Plaintext>& A,
    const std::vector<seal::Ciphertext>& B, std::vector<seal::Ciphertext>& C,
    const size_t k1, const size_t k2, const size_t k3,
    std::shared_ptr<SmallPContext> smallpcontext,
    std::shared_ptr<BigPContext> bigpcontext = nullptr);

void matrix_mult_general_ntt_fast(const std::vector<seal::Plaintext>& A,
                                  const std::vector<seal::Ciphertext>& B,
                                  std::vector<seal::Ciphertext>& C,
                                  const size_t k1, const size_t k2,
                                  const size_t k3, seal::Evaluator& evaluator,
                                  seal::BatchEncoder& batch_encoder,
                                  seal::GaloisKeys& gal_keys);

void matrix_mult_general_fast_ctct(const std::vector<seal::Ciphertext>& A,
                                   const std::vector<seal::Ciphertext>& B,
                                   std::vector<seal::Ciphertext>& C,
                                   const size_t k1, const size_t k2,
                                   const size_t k3,
                                   std::shared_ptr<BigPContext> smallpcontext,
                                   std::shared_ptr<BigPContext> bigpcontext);

/*
 * Wrapper function around matrix multiplication. Sets up the HE context
 * and generate sample matrices of the given size.
 *
 * @param k1, k2, k3: Matrix dimensions in units of YONGSOO_LIMIT i.e.,
 *                    actual matrix dimensions are
 *                    (k1*YONGSOO_LIMIT x k2*YONGSOO_LIMIT) and
 *                    (k2*YONGSOO_LIMIT x k3*YONGSOO_LIMIT)
 * @param option3:    option3 set to true will run plaintext-ciphertext
 *                    matrix multiplication. Set to false, it will run
 *                    ciphertext-ciphertext matrix multiplication
 */
void mult_wrapper(const size_t k1, const size_t k2, const size_t k3,
                  const size_t selection, const size_t logp);

void mult_wrapper_bigp(const size_t k1, const size_t k2, const size_t k3,
                       const size_t selection, const size_t logp);

/*
 * Helper function: Checks if passed number is a power of 2
 */
inline bool power_of_two(size_t d) { return (d and (!(d & (d - 1)))); }

/*
 * Helper function: Prints the name of the example in a fancy banner.
 */
inline void print_example_banner(std::string title) {
  if (!title.empty()) {
    std::size_t title_length = title.length();
    std::size_t banner_length = title_length + 2 * 10;
    std::string banner_top = "+" + std::string(banner_length - 2, '-') + "+";
    std::string banner_middle =
        "|" + std::string(9, ' ') + title + std::string(9, ' ') + "|";

    std::cout << std::endl
              << banner_top << std::endl
              << banner_middle << std::endl
              << banner_top << std::endl;
  }
}

/*
 * Helper function: Prints the `parms_id' to std::ostream.
 */
inline std::ostream& operator<<(std::ostream& stream,
                                seal::parms_id_type parms_id) {
  /*
   Save the formatting information for std::cout.
   */
  std::ios old_fmt(nullptr);
  old_fmt.copyfmt(std::cout);

  stream << std::hex << std::setfill('0') << std::setw(16) << parms_id[0] << " "
         << std::setw(16) << parms_id[1] << " " << std::setw(16) << parms_id[2]
         << " " << std::setw(16) << parms_id[3] << " ";

  /*
   Restore the old std::cout formatting.
   */
  std::cout.copyfmt(old_fmt);

  return stream;
}

/*
 * Helper function: C = A x B where A and B are square
 */

template <class T>
inline void mul(std::vector<T> A, std::vector<T> B, std::vector<T>& C) {
  long dim = sqrt(A.size());

  for (long i = 0; i < dim; i++) {
    for (long j = 0; j < dim; j++) {
      T acc = T(0);
      for (long k = 0; k < dim; k++) {
        T tmp = A[i * dim + k] * B[k * dim + j];
        acc += tmp;
      }
      C[i * dim + j] = acc;
    }
  }
}

template <class T>
inline void mul(std::vector<T> A, std::vector<T>& B, T k) {
  for (long i = 0; i < A.size(); i++) {
    B[i] = A[i] * k;
  }
}

template <class T>
inline bool check_equal(std::vector<T>& a, std::vector<T>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

// Wrapper for ct-ct mul
inline void ct_mul(const seal::Ciphertext& a, const seal::Ciphertext& b,
                   seal::Ciphertext& c, seal::Evaluator& evaluator,
                   std::shared_ptr<BigPContext> bigpcontext = nullptr) {
  if (bigpcontext != nullptr) {
    ct_ct_mul(a, b, c, bigpcontext);
  } else {
    evaluator.multiply(a, b, c);
  }
}

// Wrapper for ct-pt mul
inline void ct_pt_mul(const seal::Ciphertext& a, const seal::Plaintext& b,
                      seal::Ciphertext& c, seal::Evaluator& evaluator,
                      std::shared_ptr<BigPContext> bigpcontext = nullptr) {
  if (bigpcontext != nullptr) {
    c = a;
    pt_ct_mul(c, b, bigpcontext);
  } else {
    evaluator.multiply_plain(a, b, c);
  }
}
