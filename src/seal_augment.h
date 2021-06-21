#pragma once

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"

#include <assert.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include <seal/seal.h>

#define GALKEY true
#define FASTBLOCK true
#define HOISTING true
#define YONGSOO_LIMIT (128)
#define TIME_UNIT milliseconds
#define STRING(X) #X
#define PRINT_STRING(X) " " STRING(X)

void test_ntt(int logp, int logn);

/*
 * Helper function: Prints the parameters in a SEALContext.
 */
inline void print_parameters(std::shared_ptr<seal::SEALContext> context) {
  // Verify parameters
  if (!context) {
    throw std::invalid_argument("context is not set");
  }
  auto &context_data = *context->key_context_data();

  /*
   Which scheme are we using?
   */
  std::string scheme_name;
  switch (context_data.parms().scheme()) {
    case seal::scheme_type::BFV:
      scheme_name = "BFV";
      break;
    case seal::scheme_type::CKKS:
      scheme_name = "CKKS";
      break;
    default:
      throw std::invalid_argument("unsupported scheme");
  }
  std::cout << "/" << std::endl;
  std::cout << "| Encryption parameters :" << std::endl;
  std::cout << "|   scheme: " << scheme_name << std::endl;
  std::cout << "|   poly_modulus_degree: "
            << context_data.parms().poly_modulus_degree() << std::endl;

  /*
   Print the size of the true (product) coefficient modulus.
   */
  std::cout << "|   coeff_modulus size: ";
  std::cout << context_data.total_coeff_modulus_bit_count() << " (";
  auto coeff_modulus = context_data.parms().coeff_modulus();
  std::size_t coeff_mod_count = coeff_modulus.size();
  for (std::size_t i = 0; i < coeff_mod_count - 1; i++) {
    std::cout << coeff_modulus[i].bit_count() << " + ";
  }
  std::cout << coeff_modulus.back().bit_count();
  std::cout << ") bits" << std::endl;

  /*
   For the BFV scheme print the plain_modulus parameter.
   */
  if (context_data.parms().scheme() == seal::scheme_type::BFV) {
    std::cout << "|   plain_modulus: "
              << context_data.parms().plain_modulus().value() << std::endl;
  }

  std::cout << "\\" << std::endl;
}

struct SEALStuff {
  std::shared_ptr<seal::SEALContext> context;
  std::shared_ptr<seal::EncryptionParameters> parms;

  std::shared_ptr<seal::Evaluator> evaluator;
};

struct SEALStuffClient {
  std::unique_ptr<seal::KeyGenerator> key_generator;
  seal::PublicKey public_key;
  seal::SecretKey secret_key;
  std::unique_ptr<seal::Encryptor> encryptor;
  std::unique_ptr<seal::Decryptor> decryptor;
  seal::GaloisKeys gal_keys;
  seal::RelinKeys relin_keys;
};

struct BigPContext {
  std::unique_ptr<SEALStuff> seal_stuff;
  std::unique_ptr<SEALStuffClient> seal_stuff_client;
  int32_t n;
  NTL::ZZ_p n_inv;
  NTL::ZZ p;
  std::vector<NTL::ZZ> qs;
  NTL::ZZ q;
  NTL::ZZ_p xi;
  NTL::ZZ_p xi_inv;
  NTL::ZZ_p omega;
  std::vector<NTL::ZZ_p> omega_pow;
  std::vector<int> index_map;

  std::vector<uint64_t> mods1, mods2;
};

struct SmallPContext {
  std::unique_ptr<SEALStuff> seal_stuff;
  std::unique_ptr<SEALStuffClient> seal_stuff_client;
  std::unique_ptr<seal::BatchEncoder> batch_encoder;
  // int32_t n;
  // NTL::ZZ p;
  // std::vector<NTL::ZZ> qs;
  // NTL::ZZ q;
  // NTL::ZZ_p xi;
  // NTL::ZZ_p xi_inv;
  // NTL::ZZ_p omega;
  // std::vector<int> index_map;
  // std::vector<uint64_t> mods1, mods2;
};

std::shared_ptr<BigPContext> setup_params(size_t logp, size_t logn,
                                          size_t num_q, size_t logqi,
                                          int selection);

std::shared_ptr<SmallPContext> setup_params(size_t logp, size_t logn,
                                            int selection);

seal::Plaintext big_encode(const std::vector<NTL::ZZ_p> &plaintext,
                           std::shared_ptr<BigPContext> context);

// plaintext is something like m mod qi for each i.
std::vector<NTL::ZZ_p> big_decode(const seal::Plaintext &plaintext,
                                  std::shared_ptr<BigPContext> context);

NTL::ZZ from_uint64_t(const uint64_t &x);

uint64_t to_uint64_t(const NTL::ZZ &x);

NTL::ZZ_p to_zzp(const NTL::ZZ &x);

seal::Ciphertext big_encrypt(const std::vector<NTL::ZZ_p> &plaintext,
                             std::shared_ptr<BigPContext> context);

std::vector<NTL::ZZ_p> big_decrypt(const seal::Ciphertext &ciphertext,
                                   std::shared_ptr<BigPContext> context);

void rotate(seal::Ciphertext &ciphertext, int32_t step,
            std::shared_ptr<BigPContext> context);

void pt_ct_mul(seal::Ciphertext &ct, const seal::Plaintext &pt,
               std::shared_ptr<BigPContext> context);

void ct_ct_mul(const seal::Ciphertext &a, const seal::Ciphertext &b,
               seal::Ciphertext &c, std::shared_ptr<BigPContext> context);

void test_encode_decode();
void test_encrypt_decrypt();
void test_multiply();

std::vector<NTL::ZZ_p> ntt_inv_fast(const std::vector<NTL::ZZ_p> &a,
                                    int32_t log_n, const NTL::ZZ_p &omega);

std::vector<NTL::ZZ_p> ntt_fast(const std::vector<NTL::ZZ_p> &a, int32_t log_n,
                                const NTL::ZZ_p &omega);

void ntt_superfast(std::vector<NTL::ZZ_p> &a, int32_t log_n,
                   const std::vector<NTL::ZZ_p> &omega_pow, int32_t mul);

void ntt_inv_superfast(std::vector<NTL::ZZ_p> &a, int32_t log_n,
                       const std::vector<NTL::ZZ_p> &omega_pow,
                       const NTL::ZZ_p &n_inv);

void test_rotate(std::shared_ptr<BigPContext> context);


void test_linear_monomials(size_t u , size_t v);