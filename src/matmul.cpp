/*******************************************************
 * 2019 Sameer Wagh <snwagh@gmail.com>
 *
 * Ciphertext -ciphertext matrix multiplication using SEAL
 * and BFV implemented. This assumes columnwise encoding of
 * matrices. Algorithm used is CCS18 paper with
 * several optimization techniques.
 *
 * This file is part of ponytail project.
 *******************************************************/

#include "matmul.h"
#include <NTL/BasicThreadPool.h>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;
using namespace seal;
using namespace NTL;

size_t SUB_OPTIMAL;
size_t NUM_THREADS;
const long memoryscale = (1 << 20);   //! 2^20: linux, 2^30: mac

#define TIMING true
#define parallel false
/*
 * Generate \tilde{A} and \tilde{B} of size slot_count each.
 *
 * Inputs:  Matrices A, B of size d x d each with some data.
 * Outputs: Data contains 2 or 4 matrices (depending on SUB_OPTIMAL factor)
 *          If 2 then both top and bottom row contain A, B as YONGSOO_LIMIT
 *          size matrices which correspond to the \tilde{A} and \tilde{B}
 *          assuming A and B were padded with zeros. If 4, then the two
 *          matrices in the top and bottom rows are interleaved.
 */

template <class T>
void generate_tilde(vector<T> &A, vector<T> &B, vector<T> &data_1,
                    vector<T> &data_2, size_t d, size_t slot_count) {
  assert((A.size() == d * d and B.size() == d * d) && "Size mismatch A, B");
  assert((data_1.size() == slot_count and data_2.size() == slot_count) &&
         "Size mismatch data_1, data_2");

  vector<T> temp_1(YONGSOO_LIMIT * YONGSOO_LIMIT);
  size_t temp_size = SUB_OPTIMAL * YONGSOO_LIMIT * YONGSOO_LIMIT;

  /*****************************************************************
   Fit small matrix into a standard YONGSOO_LIMIT size matrix and then
   use do the initial transformation on the YONGSOO_LIMIT matrix
   *****************************************************************/
  for (size_t i = 0; i < d; ++i)
    for (size_t j = 0; j < d; ++j) temp_1[i * YONGSOO_LIMIT + j] = A[i * d + j];

  if (SUB_OPTIMAL == 1) {
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i)
      for (size_t j = 0; j < YONGSOO_LIMIT; ++j) {
        data_1[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
        data_1[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
      }
  } else {
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i)
      for (size_t j = 0; j < YONGSOO_LIMIT; ++j) {
        data_1[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
        data_1[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j) + 1] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
        data_1[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
        data_1[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j) + 1] =
            temp_1[i * YONGSOO_LIMIT + ((j + i) % YONGSOO_LIMIT)];
      }
  }

  vector<T> temp_2(YONGSOO_LIMIT * YONGSOO_LIMIT);
  for (size_t i = 0; i < d; ++i)
    for (size_t j = 0; j < d; ++j) temp_2[i * YONGSOO_LIMIT + j] = B[i * d + j];

  if (SUB_OPTIMAL == 1) {
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i)
      for (size_t j = 0; j < YONGSOO_LIMIT; ++j) {
        data_2[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
        data_2[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
      }
  } else {
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i)
      for (size_t j = 0; j < YONGSOO_LIMIT; ++j) {
        data_2[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
        data_2[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j) + 1] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
        data_2[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
        data_2[temp_size + SUB_OPTIMAL * (i * YONGSOO_LIMIT + j) + 1] =
            temp_2[((i + j) % YONGSOO_LIMIT) * YONGSOO_LIMIT + j];
      }
  }

#if (VERBOSE)
  cout << endl;
  for (size_t rows = 0; rows < YONGSOO_LIMIT; ++rows) {
    for (size_t cols = 0; cols < YONGSOO_LIMIT; ++cols) {
      T val = data_1[SUB_OPTIMAL * (rows * YONGSOO_LIMIT + cols)];
      cout << "(" << (T)val / YONGSOO_LIMIT << "," << val % YONGSOO_LIMIT << ")"
           << "\t";
    }
    cout << endl;
  }
#endif
}

template <class T>
void decode_matrix(vector<T> &A, const vector<T> B, size_t dim) {
  for (size_t rows = 0; rows < dim; ++rows) {
    for (size_t cols = 0; cols < dim; ++cols) {
      A[rows * dim + cols] = B[SUB_OPTIMAL * (rows * YONGSOO_LIMIT + cols)];
    }
  }
}

/*
 * Generate right shifting masks for matrix A.
 * These masks are of size YONGSOO_LIMIT each.
 */
void generate_masks(vector<Plaintext> &right_mask,
                    BatchEncoder &batch_encoder) {
  size_t slot_count = batch_encoder.slot_count();

  vector<vector<uint64_t>> masks_r(YONGSOO_LIMIT,
                                   vector<uint64_t>(slot_count, 0));

  assert((masks_r[1][0] == 0) && "Masks must be zeroed out before passing");
  size_t temp_size = 2 * SUB_OPTIMAL * YONGSOO_LIMIT * YONGSOO_LIMIT;

  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    for (size_t j = 0; j < temp_size; ++j) {
      if ((j % (SUB_OPTIMAL * YONGSOO_LIMIT)) < (SUB_OPTIMAL * i))
        masks_r[i][j] = 1;
    }
    batch_encoder.encode(masks_r[i], right_mask[i]);
  }
}

/*
 * Generate left and right shifting masks for matrix A.
 * These masks are of size YONGSOO_LIMIT each.
 */
template <class T>
void generate_masks_hoisting(vector<vector<T>> &right_mask,
                             vector<vector<T>> &left_mask, size_t slot_count) {
  right_mask.resize(YONGSOO_LIMIT);
  left_mask.resize(YONGSOO_LIMIT);
  for (int i = 0; i < YONGSOO_LIMIT; i++) {
    right_mask[i].resize(slot_count);
    left_mask[i].resize(slot_count);
  }

  assert((right_mask[1][0] == 0 and left_mask[1][0] == 0) &&
         "Masks must be zeroed out before passing");
  size_t temp_size = 2 * SUB_OPTIMAL * YONGSOO_LIMIT * YONGSOO_LIMIT;

  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    for (size_t j = 0; j < temp_size; ++j) {
      size_t j1 = (j % (SUB_OPTIMAL * YONGSOO_LIMIT));
      size_t q1 = floor(j / (SUB_OPTIMAL * YONGSOO_LIMIT)) + 1;
      size_t q2 = q1 * (SUB_OPTIMAL * YONGSOO_LIMIT) - 1;
      size_t j2 = q2 - j1;  // shifted position

      if (j1 < (SUB_OPTIMAL * i))
        right_mask[i][j2] = T(1);
      else
        left_mask[i][j2] = T(1);
    }
  }
}

// The easiest method is to just apply rotation all together
// like evaluator.rotate_rows_many(ct, steps, gal_keys, res);
// If we can use multi-cores, then divide into several blocks and rotate for
// each block
void rotate_rows_many_divide(Ciphertext ct, vector<int> steps,
                             Evaluator &evaluator, GaloisKeys &gal_keys,
                             vector<Ciphertext> &res) {
  res[0] = ct;
  if (NUM_THREADS == 1) {
    evaluator.rotate_rows_many(res[0], steps, gal_keys, res);
  }
  else {
    size_t total_len = steps.size();  // the total number of needed rotations
    size_t block_len = (size_t)ceil((double)total_len / NUM_THREADS);
    size_t final_block_len = (total_len % block_len);
    if (final_block_len == 0) final_block_len = block_len;

    NTL_EXEC_RANGE(NUM_THREADS, first, last);
    for (int k = first; k < last; ++k) {
      size_t nstart = k * block_len;
      size_t sub_len = (k == (NUM_THREADS - 1) ? final_block_len : block_len);

      vector<int> new_steps;
      for (size_t i = nstart; i < nstart + sub_len; ++i) {
        new_steps.push_back(steps[i]);
      }

      vector<Ciphertext> ct_temp(sub_len + 1);
      evaluator.rotate_rows_many(res[0], new_steps, gal_keys, ct_temp);

      for (size_t i = 1; i < sub_len + 1; ++i) {
        res[i + nstart] = ct_temp[i];
      }
    }
    NTL_EXEC_RANGE_END
  }
}

/*
 * Generate ciphertext components from the initial ciphertexts.
 * A_cts[0] should contain original tilde ciphertexts i.e., rotations.
 */
void generate_components_Acts(vector<Ciphertext> &A_cts,
                              vector<Plaintext> &right_mask,
                              Evaluator &evaluator, GaloisKeys &gal_keys) {
  assert(A_cts.size() == YONGSOO_LIMIT && "Size mismatch for generating cts");

  NTL_EXEC_RANGE(YONGSOO_LIMIT - 1, first, last);
  for (int i = first; i < last; ++i) {
    size_t i1 = i + 1;
    Ciphertext A_temp_r = A_cts[0];
    Ciphertext A_temp_l = A_cts[0];
    evaluator.multiply_plain_inplace(A_temp_r, right_mask[i1]);
    evaluator.sub_inplace(A_temp_l, A_temp_r);
    evaluator.rotate_rows_inplace(A_temp_r, SUB_OPTIMAL * (-YONGSOO_LIMIT + i1),
                                  gal_keys);
    evaluator.rotate_rows_inplace(A_temp_l, SUB_OPTIMAL * i1, gal_keys);
    evaluator.add(A_temp_r, A_temp_l, A_cts[i1]);
  }
  NTL_EXEC_RANGE_END
}

void generate_components_hoisting_Acts(
    vector<Ciphertext> &A_cts, vector<Plaintext> &right_mask,
    vector<Plaintext> &left_mask, Evaluator &evaluator, GaloisKeys &gal_keys,
    std::shared_ptr<BigPContext> bigpcontext = nullptr) {
  assert(A_cts.size() == YONGSOO_LIMIT && "Size mismatch for generating cts");

  vector<int> steps;
  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    steps.push_back(SUB_OPTIMAL * (-YONGSOO_LIMIT + i));
    steps.push_back(SUB_OPTIMAL * i);
  }

  vector<Ciphertext> A_temp(2 * YONGSOO_LIMIT - 1);
  A_temp[0] = A_cts[0];
  rotate_rows_many_divide(
      A_temp[0], steps, evaluator, gal_keys,
      A_temp);  // this part has been changed for efficient multi-threading mode

  NTL_EXEC_RANGE(YONGSOO_LIMIT - 1, first, last);
  for (int i = first; i < last; ++i) {
    size_t i1 = i + 1;
    Ciphertext A_temp_r = A_temp[2 * i1 - 1];
    Ciphertext A_temp_l = A_temp[2 * i1];

    if (bigpcontext != nullptr) {
      pt_ct_mul(A_temp_r, right_mask[i1], bigpcontext);
      pt_ct_mul(A_temp_l, left_mask[i1], bigpcontext);
    } else {  // small p version
      evaluator.multiply_plain_inplace(A_temp_r, right_mask[i1]);
      evaluator.multiply_plain_inplace(A_temp_l, left_mask[i1]);
    }
    evaluator.add(A_temp_r, A_temp_l, A_cts[i1]);
  }
  NTL_EXEC_RANGE_END
}

template <class T>
void prepare_vec_helper(vector<T> &input, vector<T> &output, size_t k) {
  size_t ROW_TWO = SUB_OPTIMAL * YONGSOO_LIMIT * YONGSOO_LIMIT;
  for (size_t i = 0; i < YONGSOO_LIMIT;
       ++i)  // SUB_OPTIMAL * YONGSOO_LIMIT; ++row)
  {
    for (size_t j = 0; j < YONGSOO_LIMIT - k; ++j) {
      if (SUB_OPTIMAL == 1) {
        T val = input[i * YONGSOO_LIMIT + (j + k)];
        output[i * YONGSOO_LIMIT + j] = val;
        output[ROW_TWO + i * YONGSOO_LIMIT + j] = val;
      } else {
        T val = input[i * YONGSOO_LIMIT * SUB_OPTIMAL + (j + k) * SUB_OPTIMAL];
        output[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL] = val;
        output[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL + 1] = val;
        output[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL] =
            val;
        output[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL +
               1] = val;
      }
    }
    for (size_t j = 0; j < k; ++j) {
      if (SUB_OPTIMAL == 1) {
        size_t j1 = j + YONGSOO_LIMIT - k;
        T val = input[i * YONGSOO_LIMIT + j];
        output[i * YONGSOO_LIMIT + j1] = val;
        output[ROW_TWO + i * YONGSOO_LIMIT + j1] = val;
      } else {
        size_t j1 = j + YONGSOO_LIMIT - k;
        T val = input[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL];
        output[i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL] = val;
        output[i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL + 1] = val;
        output[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL] =
            val;
        output[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL +
               1] = val;
      }
    }
  }
}

/*
 * Generate plaintext components from the initial plaintext.
 * A_cts[0] contain original tilde
 * (i.e., rotations) of plaintexts.
 */

void generate_components_Apts(vector<Plaintext> &A_pts,
                              BatchEncoder &batch_encoder) {
  assert(A_pts.size() == YONGSOO_LIMIT && "Size mismatch for generating cts");
  size_t slot_count = batch_encoder.slot_count();

  vector<uint64_t> data_A(slot_count), temp(slot_count);
  batch_encoder.decode(A_pts[0], data_A);

  for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
    prepare_vec_helper<uint64_t>(data_A, temp, k);
    batch_encoder.encode(temp, A_pts[k]);
  }
}
// kind of duplicate but it also works for large p.
// The above two functions could be merged into one.
void generate_components_Apts(vector<Plaintext> &A_pts,
                              std::shared_ptr<SmallPContext> smallpcontext,
                              std::shared_ptr<BigPContext> bigpcontext) {
  assert(A_pts.size() == YONGSOO_LIMIT && "Size mismatch for generating cts");

  size_t slot_count;

  if (smallpcontext != nullptr) {
    slot_count = smallpcontext->batch_encoder->slot_count();
  } else {
    slot_count = bigpcontext->seal_stuff->parms->poly_modulus_degree();  //
  }
  // vector<uint64_t> data_A(slot_count), temp(slot_count);

  // Need to take care of this step.
  if (smallpcontext != nullptr) {
    vector<uint64_t> data_A(slot_count), temp(slot_count);
    smallpcontext->batch_encoder->decode(A_pts[0], data_A);
    // need to refactor the code
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
      prepare_vec_helper<uint64_t>(data_A, temp, k);
      smallpcontext->batch_encoder->encode(temp, A_pts[k]);
    }
  } else {
    vector<NTL::ZZ_p> data_A(slot_count), temp(slot_count);
    data_A = big_decode(A_pts[0], bigpcontext);  // fixme implement

    /*--------------------------------------------*/
    // PARALLELIZATION !!!!!!!!!!!!!!!!!!!!
    /*--------------------------------------------*/
#if parallel
    NTL_EXEC_RANGE(YONGSOO_LIMIT - 1, first, last);
    for(size_t k = first; k < last; k++){
      prepare_vec_helper<NTL::ZZ_p>(data_A, temp, k + 1);
      A_pts[k + 1] = big_encode(temp, bigpcontext);
    }
    NTL_EXEC_RANGE_END
#else
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
      prepare_vec_helper<NTL::ZZ_p>(data_A, temp, k); // uncomment for correctness.
      A_pts[k] = big_encode(temp, bigpcontext);
    }
#endif
  }
}

/*
 * Generate all the ciphertext components from the initial ciphertexts.
 * A_cts[i][0] should contain original tilde ciphertexts i.e., rotations.
 */
void generate_all_components_Acts(vector<vector<Ciphertext>> &A_cts,
                                  const vector<Ciphertext> &A,
                                  vector<Plaintext> &right_mask,
                                  Evaluator &evaluator, GaloisKeys &gal_keys) {
  for (size_t i = 0; i < A.size(); ++i) {
    vector<Ciphertext> A_tmp(YONGSOO_LIMIT);
    A_tmp[0] = A[i];
    generate_components_Acts(A_tmp, right_mask, evaluator, gal_keys);
    A_cts[i] = A_tmp;
  }
}

void generate_all_components_hoisting_Acts(
    vector<vector<Ciphertext>> &A_cts, const vector<Ciphertext> &A,
    vector<Plaintext> &right_mask, vector<Plaintext> &left_mask,
    Evaluator &evaluator, GaloisKeys &gal_keys,
    shared_ptr<BigPContext> bigpcontext = nullptr) {
 
  for (size_t i = 0; i < A.size(); ++i) {
    vector<Ciphertext> A_tmp(YONGSOO_LIMIT);
    A_tmp[0] = A[i];
    generate_components_hoisting_Acts(A_tmp, right_mask, left_mask, evaluator,
                                      gal_keys, bigpcontext);
    A_cts[i] = A_tmp;
  }
}

/*
 * Generate all the plaintext components from the initial plaintexts.
 * A_cts[i][0] should contain original tilde i.e., rotations.
 */
void generate_all_components_Apts(vector<vector<Plaintext>> &A_pts,
                                  const vector<Plaintext> &A,
                                  std::shared_ptr<SmallPContext> smallpcontext,
                                  std::shared_ptr<BigPContext> bigpcontext)
// BatchEncoder &batch_encoder)
{
  for (size_t i = 0; i < A.size(); ++i) {
    vector<Plaintext> A_tmp(YONGSOO_LIMIT);
    A_tmp[0] = A[i];
    generate_components_Apts(A_tmp, smallpcontext, bigpcontext);
    A_pts[i] = A_tmp;
  }
}

/*
 * Generate ciphertext components from the initial ciphertexts.
 * B_cts[0] should contain original tilde ciphertexts i.e., rotations.
 */
void generate_components_Bcts(vector<Ciphertext> &B_cts, Evaluator &evaluator,
                              GaloisKeys &gal_keys) {
  // Rotations in B are free due to repetitions
  NTL_EXEC_RANGE(YONGSOO_LIMIT - 1, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.rotate_rows(B_cts[0], (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT,
                          gal_keys, B_cts[i + 1]);
  }
  NTL_EXEC_RANGE_END
}

// This should be okay for any p. Small or large.
void generate_components_hoisting_Bcts(vector<Ciphertext> &B_cts,
                                       Evaluator &evaluator,
                                       GaloisKeys &gal_keys) {
  // Rotations in B are free due to repetitions
  vector<int> steps(YONGSOO_LIMIT - 1, 0);
  for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
    steps[i] = (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT;
  }

  rotate_rows_many_divide(
      B_cts[0], steps, evaluator, gal_keys,
      B_cts);  // this part has been changed for efficient multi-threading mode
}

/*
 * Generate all the ciphertext components from the initial ciphertexts.
 * B_cts[i][0] should contain original tilde ciphertexts i.e., rotations.
 */
void generate_all_components_Bcts(vector<vector<Ciphertext>> &B_cts,
                                  const vector<Ciphertext> &B,
                                  Evaluator &evaluator, GaloisKeys &gal_keys) {
  for (size_t i = 0; i < B.size(); ++i) {
    vector<Ciphertext> B_tmp(YONGSOO_LIMIT);
    B_tmp[0] = B[i];
#if (HOISTING)
    generate_components_hoisting_Bcts(B_tmp, evaluator, gal_keys);
#else
    generate_components_Bcts(B_tmp, evaluator, gal_keys);
#endif
    B_cts[i] = B_tmp;
  }
  // NTL_EXEC_RANGE_END
}

/*
 * Take inputs as all the ciphertext components from the initial ciphertexts of
 * A (whole matrix) and B (one column) perform blockwise matrix multiplication
 * as A * B_col
 */

void matrix_mult_unit_colwise(const vector<vector<Ciphertext>> &A_cts,
                              const vector<vector<Ciphertext>> &B_col_cts,
                              vector<Ciphertext> &C, Evaluator &evaluator,
                              RelinKeys &relin_keys,
                              std::shared_ptr<BigPContext> bigpcontext) {
  size_t k1 = C.size();
  size_t k2 = A_cts.size() / k1;
  assert(k2 == B_col_cts.size());

  // NTL_EXEC_RANGE(k1, first, last);
  // for (int i = first; i < last; ++i)
  for (size_t i = 0; i < k1; i++) {
    Ciphertext res;
    res.reserve(3);

    // new version for MT
    vector<Ciphertext> temp(YONGSOO_LIMIT);

    NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last)
    for (size_t k = first; k < last; k++) {
      ct_mul(A_cts[i * k2][k], B_col_cts[0][k], temp[k], evaluator,
             bigpcontext);
    }
    NTL_EXEC_RANGE_END

    res = temp[0];
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k)
      evaluator.add_inplace(res, temp[k]);

    // A[i * k2 + j][-] * B[j][-]
    for (size_t j = 1; j < k2; ++j) {
      NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last)
      for (size_t k = first; k < last; k++) {
        ct_mul(A_cts[i * k2 + j][k], B_col_cts[j][k], temp[k], evaluator,
               bigpcontext);
      }
      NTL_EXEC_RANGE_END

      for (size_t k = 0; k < YONGSOO_LIMIT; ++k)
        evaluator.add_inplace(res, temp[k]);
    }

    evaluator.relinearize_inplace(res, relin_keys);
    C[i] = res;
  }
}


/*
 * Take inputs as all the plaintext/ciphertext components from the initial
 * ciphertexts of A (one row) and B (whole matrix) perform blockwise matrix
 * multiplication as A_row * B
 */
void matrix_mult_unit_rowwise(const vector<vector<Plaintext>> &A_row_pts,
                              const vector<vector<Ciphertext>> &B_cts,
                              vector<Ciphertext> &C, Evaluator &evaluator,
                              std::shared_ptr<BigPContext> bigpcontext) {
  size_t k3 = C.size();
  size_t k2 = B_cts.size() / k3;
  assert(k2 == A_row_pts.size());

  for (int j = 0; j < k3; ++j) {
    Ciphertext res;
    res.reserve(3);

    // new version for MT
    vector<Ciphertext> temp(YONGSOO_LIMIT);
    NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
    for (size_t k = first; k < last; k++) {
      ct_pt_mul(B_cts[j][k], A_row_pts[0][k], temp[k], evaluator, bigpcontext);
    }
    NTL_EXEC_RANGE_END

    res = temp[0];
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k)
      evaluator.add_inplace(res, temp[k]);

    for (size_t i = 1; i < k2; ++i) {
      NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
      for (size_t k = first; k < last; k++) {
        ct_pt_mul(B_cts[i * k3 + j][k], A_row_pts[i][k], temp[k], evaluator,
                  bigpcontext);
      }
      NTL_EXEC_RANGE_END

      for (size_t k = 0; k < YONGSOO_LIMIT; ++k)
        evaluator.add_inplace(res, temp[k]);
    }

    C[j] = res;
  }
}

void matrix_mult_general_fast_ctct(const vector<Ciphertext> &A,
                                   const vector<Ciphertext> &B,
                                   vector<Ciphertext> &C, const size_t k1,
                                   const size_t k2, const size_t k3,
                                   std::shared_ptr<SmallPContext> smallpcontext,
                                   std::shared_ptr<BigPContext> bigpcontext) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

#if TIMING
  chrono::high_resolution_clock::time_point time_start, time_end;
  chrono::TIME_UNIT time_total(0);
  chrono::TIME_UNIT time_total_genA(0);
  chrono::TIME_UNIT time_total_genB(0);
  chrono::TIME_UNIT time_total_mult(0);
  time_start = chrono::high_resolution_clock::now();
  struct rusage usage;
#endif

  GaloisKeys gal_keys;
  std::shared_ptr<seal::Evaluator> evaluator;
  RelinKeys relin_keys;

  if (smallpcontext != nullptr) {
    evaluator = smallpcontext->seal_stuff->evaluator;
    gal_keys = smallpcontext->seal_stuff_client->gal_keys;
    relin_keys = smallpcontext->seal_stuff_client->relin_keys;
  } else {
    evaluator = bigpcontext->seal_stuff->evaluator;
    gal_keys = bigpcontext->seal_stuff_client->gal_keys;
    relin_keys = bigpcontext->seal_stuff_client->relin_keys;
  }
    
#if TIMING
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "> Key copy time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
       << endl;
  int ret = getrusage(RUSAGE_SELF, &usage);
  cout << "Memory Usage: \t\t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
  cout << "============================================" << endl;
  time_start = chrono::high_resolution_clock::now();
#endif
    
  /*-----------------------------------*/
  // Step 0. Generation of Masking
  /*-----------------------------------*/

#if (HOISTING)
  cout << "USing hoisting";
  vector<Plaintext> right_mask(YONGSOO_LIMIT), left_mask(YONGSOO_LIMIT);
  // generate_masks_hoisting<uint64_t>(masks_r, masks_l, batch_encoder);

  if (smallpcontext != nullptr) {
    cout << " with SMall P " << endl;
    vector<vector<uint64_t>> masks_r(YONGSOO_LIMIT), masks_l(YONGSOO_LIMIT);
    generate_masks_hoisting<uint64_t>(
        masks_r, masks_l, smallpcontext->batch_encoder->slot_count());
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i) {
      smallpcontext->batch_encoder->encode(masks_r[i], right_mask[i]);
      smallpcontext->batch_encoder->encode(masks_l[i], left_mask[i]);
    }
  } else {
    cout << " with Large P: " ;
    vector<vector<NTL::ZZ_p>> masks_r(YONGSOO_LIMIT), masks_l(YONGSOO_LIMIT);
    generate_masks_hoisting<NTL::ZZ_p>(masks_r, masks_l, bigpcontext->n);
    cout << "masking vectors gen" << endl;

    /*--------------------------------------------*/
    // PARALLELIZATION !!!!!!!!!!!!!!!!!!!!
    /*--------------------------------------------*/
#if parallel
    NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
    for (size_t i = first; i < last; ++i)
    {
      right_mask[i] = big_encode(masks_r[i], bigpcontext);
      left_mask[i] = big_encode(masks_l[i], bigpcontext);
    }
    NTL_EXEC_RANGE_END
#else
    for (size_t i = 0; i < YONGSOO_LIMIT; ++i)
    {
      right_mask[i] = big_encode(masks_r[i], bigpcontext);
      left_mask[i] = big_encode(masks_l[i], bigpcontext);
    }
#endif
  }
#else
  cout << "Not using hoisting " << endl;
  vector<Plaintext> right_mask(YONGSOO_LIMIT);
  generate_masks(right_mask, *smallpcontext->batch_encoder);
#endif

#if TIMING
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "> Mask Generation time: \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;
  ret = getrusage(RUSAGE_SELF, &usage);
  cout << "Memory Usage: \t\t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
  cout << "============================================" << endl;
  time_start = chrono::high_resolution_clock::now();
#endif

  /*-----------------------------------*/
  // Step 1. matrix-mult
  /*-----------------------------------*/

  // perform matrix-mult for each column
  if (k1 <= k3) {
    std::cout << "Got here: k1 <= k3 " << std::endl;

    /*-----------------------------------*/
    // Generate all the components of Amat
    /*-----------------------------------*/
#if (HOISTING)
    vector<vector<Ciphertext>> A_cts(k1 * k2,
                                     vector<Ciphertext>(YONGSOO_LIMIT));

    generate_all_components_hoisting_Acts(A_cts, A, right_mask, left_mask,
                                          *evaluator, gal_keys, bigpcontext);
#else
    vector<vector<Ciphertext>> A_cts(k1 * k2,
                                     vector<Ciphertext>(YONGSOO_LIMIT));
    generate_all_components_Acts(A_cts, A, right_mask, evaluator, gal_keys);
#endif
      
#if TIMING
    time_end = chrono::high_resolution_clock::now();
    time_total_genA=
        chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
    cout << "> Gen(E(A)) time: \t\t" << time_total_genA.count() << PRINT_STRING(TIME_UNIT) << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage: \t\t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "============================================" << endl;
#endif
    
    /*-----------------------------------*/
    // perform block-wise computation for each column block of matrix B
    // choose each column and store it as a vector of ciphertexts
    /*-----------------------------------*/
      
    for (size_t columns = 0; columns < k3; ++columns) {
#if TIMING
      time_start = chrono::high_resolution_clock::now();
#endif
      cout << "> START Acts(all) + Bcts(" << columns << "-col) ... " << endl;
      
      vector<Ciphertext> B_col;
      for (size_t rows = 0; rows < k2; ++rows) {
        B_col.push_back(B[rows * k3 + columns]);
      }
      vector<vector<Ciphertext>> B_col_cts(k2,
                                           vector<Ciphertext>(YONGSOO_LIMIT));
      generate_all_components_Bcts(B_col_cts, B_col, *evaluator, gal_keys);

#if TIMING
      time_end = chrono::high_resolution_clock::now();
      time_total =
          chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
      cout << "Gen(E(B" << columns << ")) time: \t\t" << time_total.count()
           << PRINT_STRING(TIME_UNIT) << endl;
      ret = getrusage(RUSAGE_SELF, &usage);
      cout<< "Memory Usage: \t\t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
      cout << "------------------------" << endl;
      time_total_genB += time_total;
      time_start = chrono::high_resolution_clock::now();
#endif
      // perform matrix computation at one time
      vector<Ciphertext> res(k1);
      matrix_mult_unit_colwise(A_cts, B_col_cts, res, *evaluator, relin_keys,
                               bigpcontext);

#if TIMING
      time_end = chrono::high_resolution_clock::now();
      time_total =
          chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
      cout << "E(A) * E(B" << columns << ") time: \t\t" << time_total.count()
           << PRINT_STRING(TIME_UNIT) << endl;
      ret = getrusage(RUSAGE_SELF, &usage);
      cout<< "Memory Usage: \t\t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
      time_total_mult += time_total;
      cout << "-----------------------------------" << endl;
#endif
      // store the computation results
      for (size_t rows = 0; rows < k1; ++rows) {
        C[k3 * rows + columns] = res[rows];
      }
    }
#if TIMING
    cout << "> total gen(B) time: \t\t" << time_total_genB.count() << PRINT_STRING(TIME_UNIT) << endl;
    cout << "> total mult time: \t\t" << time_total_mult.count() << PRINT_STRING(TIME_UNIT) << endl;
    cout << "> total (- mask) time: \t\t" << (time_total_genB +  time_total_mult + time_total_genA).count() / 2 << PRINT_STRING(TIME_UNIT)  << endl;

#endif
  }
  else {
    cout << "Got here: k1 > k3 " << endl;

    // generate all the components of Bmat
    vector<vector<Ciphertext>> B_cts(k2 * k3,
                                     vector<Ciphertext>(YONGSOO_LIMIT));
    generate_all_components_Bcts(B_cts, B, *evaluator, gal_keys);

    NTL_EXEC_RANGE(k1, first, last);
    for (size_t rows = first; rows < last; ++rows) {
      // choose each column and store it as a vector of ciphertexts
      vector<Ciphertext> A_row;
      for (size_t columns = 0; columns < k2; ++columns) {
        A_row.push_back(A[rows * k2 + columns]);
      }
#if (HOISTING)
      vector<vector<Ciphertext>> A_row_cts(k2,
                                           vector<Ciphertext>(YONGSOO_LIMIT));
      generate_all_components_hoisting_Acts(A_row_cts, A_row, right_mask,
                                            left_mask, *evaluator, gal_keys);
#else
      vector<vector<Ciphertext>> A_row_cts(k2,
                                           vector<Ciphertext>(YONGSOO_LIMIT));
      generate_all_components_Acts(A_row_cts, A_row, right_mask, evaluator,
                                   gal_keys);
#endif

      // perform matrix computation at one time
      vector<Ciphertext> res(k3);
      matrix_mult_unit_rowwise(A_row_cts, B_cts, res, *evaluator, relin_keys);

      // store the computatio results
      for (size_t columns = 0; columns < k3; ++columns) {
        C[k3 * rows + columns] = res[columns];
      }
    }
    NTL_EXEC_RANGE_END
  }
 
}

void matrix_mult_general_fast_ptct(const vector<Plaintext> &A,
                                   const vector<Ciphertext> &B,
                                   vector<Ciphertext> &C, const size_t k1,
                                   const size_t k2, const size_t k3,
                                   std::shared_ptr<SmallPContext> smallpcontext,
                                   std::shared_ptr<BigPContext> bigpcontext) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

#if TIMING
  chrono::high_resolution_clock::time_point time_start, time_end;
  chrono::TIME_UNIT time_total(0);
  chrono::microseconds time_total_genA(0);
  chrono::microseconds time_total_genB(0);

  chrono::microseconds time_total_mult(0);
  time_start = chrono::high_resolution_clock::now();
  struct rusage usage;
#endif
    
  GaloisKeys gal_keys;
  std::shared_ptr<seal::Evaluator> evaluator;
  RelinKeys relin_keys;

  if (smallpcontext != nullptr) {
    evaluator = smallpcontext->seal_stuff->evaluator;
    gal_keys = smallpcontext->seal_stuff_client->gal_keys;
    relin_keys = smallpcontext->seal_stuff_client->relin_keys;
  } else {
    evaluator = bigpcontext->seal_stuff->evaluator;
    gal_keys = bigpcontext->seal_stuff_client->gal_keys;
    relin_keys = bigpcontext->seal_stuff_client->relin_keys;
  }

#if TIMING
  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "Key time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
      << endl;
  int ret = getrusage(RUSAGE_SELF, &usage);
  cout << "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
  cout << "===================================" << endl;
  time_start = chrono::high_resolution_clock::now();
#endif

  /*-----------------------------------*/
  // generate all ciphertext components of Bmat
  /*-----------------------------------*/

  vector<vector<Ciphertext>> B_cts(k2 * k3, vector<Ciphertext>(YONGSOO_LIMIT));
  generate_all_components_Bcts(B_cts, B, *evaluator, gal_keys);

#if TIMING
  time_end = chrono::high_resolution_clock::now();
  time_total_genB = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "Gen(B) time: \t\t" << time_total_genB.count() << PRINT_STRING(TIME_UNIT) << endl;
  ret = getrusage(RUSAGE_SELF, &usage);
  cout << "Memory Usage: \t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
  cout << "===================================" << endl;
#endif

  /*-----------------------------------*/
  // perform block-wise computation for each row block of matrix A
  // choose each column and store it as a vector of ciphertexts
  /*-----------------------------------*/

  for (int i = 0; i < k1; ++i) {
      
#if TIMING
    time_start = chrono::high_resolution_clock::now();
#endif
    
    cout << "pt (" << i << "-row) + ct (all) ... ";
    vector<Plaintext> A_row;
    for (size_t j = 0; j < k2; ++j) {
      A_row.push_back(A[i * k2 + j]);
    }
    vector<vector<Plaintext>> A_row_pts(k2, vector<Plaintext>(YONGSOO_LIMIT));
    generate_all_components_Apts(A_row_pts, A_row, smallpcontext, bigpcontext);

#if TIMING
    time_end = chrono::high_resolution_clock::now();
    time_total =
        chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
    cout << "Gen(A" << i <<  ") time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
         << endl;
    time_total_genA += time_total;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Memory Usage: \t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    time_start = chrono::high_resolution_clock::now();
#endif

    // perform matrix computation at one time
    vector<Ciphertext> res(k3);
    matrix_mult_unit_rowwise(A_row_pts, B_cts, res, *evaluator, bigpcontext);

#if TIMING
    time_end = chrono::high_resolution_clock::now();
    time_total =
        chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
    cout << "A" << i << "*E(B) time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
         << endl;
    time_total_mult += time_total;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Memory Usage: \t\t" << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "------------------------" << endl;
#endif
      
    // store the computation results
    for (size_t j = 0; j < k3; ++j) {
      C[k3 * i + j] = res[j];
    }
    cout << endl;
  }
#if TIMING
  cout << "total gen(A) = " << time_total_genA.count() << PRINT_STRING(TIME_UNIT) << endl;
  cout << "total mult = " << time_total_mult.count() << PRINT_STRING(TIME_UNIT) << endl;
  cout << "Total time minus masking " << (time_total_mult + time_total_genA + time_total_genB).count() << endl;

#endif
}

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

// Hao: big prime version
void mult_wrapper_bigp(const size_t k1, const size_t k2, const size_t k3,
                       const size_t selection, const size_t logp) {
  if (selection != 1 && selection != 3) {
    cout << " NOT IMPLEMEMENTED";
    return;
  }
  size_t poly_modulus_degree = POLY_MODULUS_DEGREE;

  size_t logn = log2(poly_modulus_degree);

  size_t numq = NUM_Q;

  size_t logqi = LOG_QI;

  chrono::high_resolution_clock::time_point time_start, time_end;
  chrono::TIME_UNIT time_total(0);
  time_start = chrono::high_resolution_clock::now();

  std::shared_ptr<BigPContext> bigpcontext =
      setup_params(logp, logn, numq, logqi, selection);

  RelinKeys relin_keys = bigpcontext->seal_stuff_client->relin_keys;
  GaloisKeys gal_keys;
  vector<int> steps(3 * (YONGSOO_LIMIT - 1), 0);
  for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
      steps[i] = (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT;
      steps[i + (YONGSOO_LIMIT - 1)] = SUB_OPTIMAL * (-YONGSOO_LIMIT + i + 1);
      steps[i + (YONGSOO_LIMIT - 1) * 2] = SUB_OPTIMAL * (i + 1);
  }
  gal_keys = bigpcontext->seal_stuff_client->key_generator->galois_keys(steps);

  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << endl;
  cout << "--------------------------------------------" << endl;
  cout << "KeyGen time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
       << endl;
  cout << "--------------------------------------------" << endl;

  SUB_OPTIMAL = (POLY_MODULUS_DEGREE / (2 * YONGSOO_LIMIT * YONGSOO_LIMIT));
  assert(SUB_OPTIMAL == 1 or
         SUB_OPTIMAL == 2 && "Some issue with sub_optimality factor");

  SetNumThreads(NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  size_t d = 128;
  vector<NTL::ZZ_p> A_data(d * d), B_data(d * d), data_1(poly_modulus_degree),
      data_2(poly_modulus_degree);

  for (size_t i = 0; i < d * d; ++i) {
    A_data[i] = NTL::random_ZZ_p();
    B_data[i] = NTL::random_ZZ_p();
  }

  size_t slot_count = poly_modulus_degree;
  generate_tilde<NTL::ZZ_p>(A_data, B_data, data_1, data_2, d, slot_count);

  //  Now we need generate_tilde.
  time_start = chrono::high_resolution_clock::now();

  Ciphertext ctxt_1, ctxt_2;
  ctxt_1 = big_encrypt(data_1, bigpcontext);
  ctxt_2 = big_encrypt(data_2, bigpcontext);

  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << endl;
  cout << "--------------------------------------------" << endl;
  cout << "Enc(A0),Enc(B0) time: \t\t" << time_total.count() << PRINT_STRING(TIME_UNIT)
         << endl;
  cout << "--------------------------------------------" << endl;
    
  
  // same matrix is being repeated many times.
  vector<Ciphertext> C(k1 * k3, Ciphertext(bigpcontext->seal_stuff->context));

  shared_ptr<SmallPContext> smallpcontext = nullptr;

  time_start = chrono::high_resolution_clock::now();

  // todo: moddify this part done.
  vector<Ciphertext> B(k2 * k3, ctxt_2);

  if (selection == 1) {
    vector<Ciphertext> A(k1 * k2, ctxt_1);
    matrix_mult_general_fast_ctct(A, B, C, k1, k2, k3, smallpcontext,
                                  bigpcontext);
  } else if (selection == 3) {
    Plaintext plain = big_encode(data_1, bigpcontext);
    vector<Plaintext> A(k1 * k2, plain);
    matrix_mult_general_fast_ptct(A, B, C, k1, k2, k3, smallpcontext,
                                  bigpcontext);
  }

  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "--------------------------------------------" << endl;
  cout << "Total time (selection=" << selection << "): \t"
       << time_total.count() << PRINT_STRING(TIME_UNIT) << endl;
  cout << "--------------------------------------------" << endl;

  vector<NTL::ZZ_p> plain_3 = big_decrypt(C[0], bigpcontext);

  // vector<NTL::ZZ_p> data_3(slot_count);
  // batch_encoder.decode(plain_3, data_3);

  // decode a matrix
  vector<NTL::ZZ_p> actual(d * d), C_data(d * d);
  decode_matrix<NTL::ZZ_p>(actual, plain_3, d);

  // plain matrix
  mul<NTL::ZZ_p>(A_data, B_data, C_data);
  mul<NTL::ZZ_p>(C_data, C_data, NTL::ZZ_p(k2));

  cout << "ok? " << check_equal<NTL::ZZ_p>(C_data, actual) << endl;
}

/*-------------------------------------------------*/
// Not used Functions
/*-------------------------------------------------*/

/*
 * Generate plaintext NTT components from the initial plaintext ciphertexts
 * A_cts[0] should contain original tilde (i.e., rotations) of plaintexts.
 */
void generate_ntt_components_Apts(vector<Plaintext> &A_pts_ntt,
                                  const Ciphertext &B_cts_ntt,
                                  Evaluator &evaluator,
                                  BatchEncoder &batch_encoder) {
  assert(A_pts_ntt.size() == YONGSOO_LIMIT &&
         "Size mismatch for generating cts");
  size_t slot_count = batch_encoder.slot_count();

  vector<uint64_t> data_A(slot_count), temp(slot_count);
  batch_encoder.decode(A_pts_ntt[0], data_A);
  evaluator.transform_to_ntt_inplace(A_pts_ntt[0], B_cts_ntt.parms_id());

  size_t ROW_TWO = SUB_OPTIMAL * YONGSOO_LIMIT * YONGSOO_LIMIT;

  for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
    for (size_t i = 0; i < YONGSOO_LIMIT;
         ++i)  // SUB_OPTIMAL * YONGSOO_LIMIT; ++row)
    {
      for (size_t j = 0; j < YONGSOO_LIMIT - k; ++j) {
        if (SUB_OPTIMAL == 1) {
          uint64_t val = data_A[i * YONGSOO_LIMIT + (j + k)];
          temp[i * YONGSOO_LIMIT + j] = val;
          temp[ROW_TWO + i * YONGSOO_LIMIT + j] = val;
        } else {
          uint64_t val =
              data_A[i * YONGSOO_LIMIT * SUB_OPTIMAL + (j + k) * SUB_OPTIMAL];
          temp[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL] = val;
          temp[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL + 1] = val;
          temp[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL] =
              val;
          temp[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL +
               1] = val;
        }
      }
      for (size_t j = 0; j < k; ++j) {
        if (SUB_OPTIMAL == 1) {
          size_t j1 = j + YONGSOO_LIMIT - k;
          uint64_t val = data_A[i * YONGSOO_LIMIT + j];
          temp[i * YONGSOO_LIMIT + j1] = val;
          temp[ROW_TWO + i * YONGSOO_LIMIT + j1] = val;
        } else {
          size_t j1 = j + YONGSOO_LIMIT - k;
          uint64_t val =
              data_A[i * YONGSOO_LIMIT * SUB_OPTIMAL + j * SUB_OPTIMAL];
          temp[i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL] = val;
          temp[i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL + 1] = val;
          temp[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL] =
              val;
          temp[ROW_TWO + i * YONGSOO_LIMIT * SUB_OPTIMAL + j1 * SUB_OPTIMAL +
               1] = val;
        }
      }
    }
    batch_encoder.encode(temp, A_pts_ntt[k]);
    evaluator.transform_to_ntt_inplace(A_pts_ntt[k], B_cts_ntt.parms_id());
  }
}

/*
 * Generate all the plaintext components from the initial plaintexts (in
 * NTT-form). A_cts[i][0] should contain original tilde i.e., rotations.
 */
void generate_all_ntt_components_Apts(vector<vector<Plaintext>> &A_pts_ntt,
                                      const vector<Plaintext> &A,
                                      const Ciphertext &B_cts_ntt,
                                      Evaluator &evaluator,
                                      BatchEncoder &batch_encoder) {
  for (size_t i = 0; i < A.size(); ++i) {
    vector<Plaintext> A_tmp(YONGSOO_LIMIT);
    A_tmp[0] = A[i];
    generate_ntt_components_Apts(A_tmp, B_cts_ntt, evaluator, batch_encoder);
    A_pts_ntt[i] = A_tmp;
  }
}

/*
 * Generate ciphertext NTT components from the initial ciphertexts.
 * B_cts[0] should contain original tilde ciphertexts i.e., rotations.
 * gen(B): 1450 ms
 */
void generate_ntt_components_Bcts(vector<Ciphertext> &B_cts_ntt,
                                  Evaluator &evaluator, GaloisKeys &gal_keys) {
  // Rotations in B are free due to repetitions
  size_t coeff_count = B_cts_ntt[0].poly_modulus_degree();
  Ciphertext Btemp = B_cts_ntt[0];

  evaluator.transform_to_ntt_inplace(B_cts_ntt[0]);

  NTL_EXEC_RANGE(YONGSOO_LIMIT - 1, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.apply_galois(
        Btemp,
        seal::util::steps_to_galois_elt((i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT,
                                        coeff_count),
        gal_keys, B_cts_ntt[i + 1]);  // galois in ctxt
    evaluator.transform_to_ntt_inplace(B_cts_ntt[i + 1]);
  }
  NTL_EXEC_RANGE_END
}

/*
 * Generate ciphertext NTT components from the initial ciphertexts
 * with hoisting technique of HS18
 * B_cts[0] should contain original tilde ciphertexts i.e., rotations.
 *
 */
void generate_ntt_components_hoisting_Bcts(vector<Ciphertext> &B_cts_ntt,
                                           Evaluator &evaluator,
                                           GaloisKeys &gal_keys) {
  // Rotations in B are free due to repetitions
  vector<int> steps(YONGSOO_LIMIT - 1, 0);
  for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
    steps[i] = (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT;
  }
  rotate_rows_many_divide(B_cts_ntt[0], steps, evaluator, gal_keys,
                          B_cts_ntt);  // this part has been changed for
                                       // efficient multi-threading mode

  // Convert the genrated ciphertext in the form of NTT
  NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.transform_to_ntt_inplace(B_cts_ntt[i]);
  }
  NTL_EXEC_RANGE_END
}

/*
 * Generate all the ciphertext components from the initial ciphertexts (in
 * NTT-form). B_cts[i][0] should contain original tilde ciphertexts i.e.,
 * rotations.
 */
void generate_all_ntt_components_Bcts(vector<vector<Ciphertext>> &B_cts,
                                      const vector<Ciphertext> &B,
                                      Evaluator &evaluator,
                                      GaloisKeys &gal_keys) {
  // NTL_EXEC_RANGE(B.size(), first, last);
  // for (long i = first; i < last; ++i)
  for (long i = 0; i < B.size(); ++i) {
    vector<Ciphertext> B_tmp(YONGSOO_LIMIT);
    B_tmp[0] = B[i];
#if (HOISTING)
    generate_ntt_components_hoisting_Bcts(B_tmp, evaluator, gal_keys);
#else
    generate_ntt_components_Bcts(B_tmp, evaluator, gal_keys);
#endif
    B_cts[i] = B_tmp;
  }
  // NTL_EXEC_RANGE_END
}

/*
 * Single ciphertext matrix multiplication of CCS18
 *
 * @param A:             Ciphertext A representing a square matrix of dimension
 * YONGSOO_LIMIT. Matrix A is values of \tilde{A} i.e., the output of the first
 * transformation completed.
 * @param B:             Ciphertext B representing the second square matrix of
 * same dimension. Note that B has to be repeated ciphertext.
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 * @param keygen:        HE KeyGenerator
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix
 */
void matrix_mult_unit(const Ciphertext &A, const Ciphertext &B, Ciphertext &C,
                      Evaluator &evaluator, BatchEncoder &batch_encoder,
                      GaloisKeys &gal_keys, RelinKeys &relin_keys) {
  size_t slot_count = batch_encoder.slot_count();
  assert(slot_count >= 2 * YONGSOO_LIMIT * YONGSOO_LIMIT &&
         "Yongsoo limit exceded");

  /*****************************************************************
   Generate ciphertexts for all A, B matrices using the left shift and
   right shift masks. Initialize the 0th components to encryptions of
   \tilde{A} and \tilde{B}
   *****************************************************************/
  vector<Ciphertext> A_cts(YONGSOO_LIMIT);
  A_cts[0] = A;

#if (HOISTING)
  vector<vector<uint64_t>> masks_r(YONGSOO_LIMIT), masks_l(YONGSOO_LIMIT);
  vector<Plaintext> right_mask(YONGSOO_LIMIT), left_mask(YONGSOO_LIMIT);
  generate_masks_hoisting<uint64_t>(masks_r, masks_l, slot_count);
  for (size_t i = 0; i < YONGSOO_LIMIT; i++) {
    batch_encoder.encode(masks_r[i], right_mask[i]);
    batch_encoder.encode(masks_l[i], left_mask[i]);
  }
  // } else{
  //   right_mask[i] =  big_encode(masks_r[i]);
  //   left_mask[i] = big_encode(masks_l[i]);
  // }
  generate_components_hoisting_Acts(A_cts, right_mask, left_mask, evaluator,
                                    gal_keys);
#else
  vector<Plaintext> right_mask(YONGSOO_LIMIT);
  generate_masks(right_mask, batch_encoder, slot_count);
  generate_components_Acts(A_cts, right_mask, evaluator, gal_keys);
#endif

  vector<Ciphertext> B_cts(YONGSOO_LIMIT);
  B_cts[0] = B;
#if (HOISTING)
  generate_components_hoisting_Bcts(B_cts, evaluator, gal_keys);
#else
  generate_components_Bcts(B_cts, evaluator, gal_keys);
#endif

  /*****************************************************************
   Complete the final summation for Yongsoo multiplication as term-by-term
   products of all the ciphertexts.
   *****************************************************************/

  vector<Ciphertext> temp_product(YONGSOO_LIMIT);
  NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.multiply(A_cts[i], B_cts[i], temp_product[i]);
  }
  NTL_EXEC_RANGE_END

  C = temp_product[0];
  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    evaluator.add_inplace(C, temp_product[i]);
  }

  evaluator.relinearize_inplace(C, relin_keys);
}

/*
 * Single ciphertext matrix multiplication of CCS18
 *
 * @param A:             Plaintext A representing a square matrix of dimension
 * YONGSOO_LIMIT. Matrix A is values of \tilde{A} i.e., the output of the first
 * transformation completed.
 * @param B:             Ciphertext B representing the second square matrix of
 * same dimension. Note that B has to be repeated ciphertext.
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 * @param keygen:        HE KeyGenerator
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix
 */
void matrix_mult_unit(const Plaintext &A, const Ciphertext &B, Ciphertext &C,
                      Evaluator &evaluator, BatchEncoder &batch_encoder,
                      GaloisKeys &gal_keys) {
  assert(batch_encoder.slot_count() >= 2 * YONGSOO_LIMIT * YONGSOO_LIMIT &&
         "Yongsoo limit exceded");

  /*****************************************************************
   Generate all components for matrices A, B using the left shift and
   right shift masks. Initialize the 0th components to encryptions of
   \tilde{A} and \tilde{B}
   *****************************************************************/
  vector<Plaintext> A_pts(YONGSOO_LIMIT);
  vector<Ciphertext> B_cts(YONGSOO_LIMIT);
  A_pts[0] = A;
  B_cts[0] = B;

  generate_components_Apts(A_pts, batch_encoder);

#if (HOISTING)
  generate_components_hoisting_Bcts(B_cts, evaluator, gal_keys);
#else
  generate_components_Bcts(B_cts, evaluator, gal_keys);
#endif

  /*****************************************************************
   Complete the final summation for Yongsoo multiplication as term-by-term
   products of all the ciphertexts.
   *****************************************************************/

  vector<Ciphertext> temp_product(YONGSOO_LIMIT);
  NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.multiply_plain(B_cts[i], A_pts[i], temp_product[i]);
  }
  NTL_EXEC_RANGE_END

  C = temp_product[0];
  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    evaluator.add_inplace(C, temp_product[i]);
  }
}

/*
 * Single ciphertext matrix multiplication but with plaintext in NTT form
 *
 * @param A:             vector of Plaintexts A representing a the YONGSOO_LIMIT
 * number of rotated Plaintexts given in NTT format.
 * @param B:             Ciphertext B representing the second square matrix of
 * same dimension. Note that B has to be repeated ciphertext.
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 * @param keygen:        HE KeyGenerator
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix
 */
void matrix_mult_ntt_unit(const Plaintext &A_ntt, const Ciphertext &B,
                          Ciphertext &C, Evaluator &evaluator,
                          BatchEncoder &batch_encoder, GaloisKeys &gal_keys) {
#if (VERBOSE)
  size_t slot_count = batch_encoder.slot_count();
  assert(slot_count >= 2 * YONGSOO_LIMIT * YONGSOO_LIMIT &&
         "Yongsoo limit exceded");
#endif
  /*****************************************************************
   Generate components in NTT space. For a cyclic row rotation to the left
   by i steps, the Galois element must be an odd integer (used with
   batching here and is equal to 3^i % M) in the interval [1, M-1],
   where M = 2*N, and N = degree(poly_modulus).
   *****************************************************************/
  vector<Plaintext> A_pts_ntt(YONGSOO_LIMIT);
  vector<Ciphertext> B_cts_ntt(YONGSOO_LIMIT);

  A_pts_ntt[0] = A_ntt;
  B_cts_ntt[0] = B;
  generate_ntt_components_Apts(A_pts_ntt, B, evaluator, batch_encoder);

#if (HOISTING)
  generate_ntt_components_hoisting_Bcts(B_cts_ntt, evaluator, gal_keys);
#else
  generate_ntt_components_Bcts(B_cts_ntt, evaluator, gal_keys);
#endif

  /*****************************************************************
   Complete the final summation for Yongsoo multiplication as term-by-term
   products of all the ciphertexts in NTT space. Final result in non NTT space.
   *****************************************************************/

  vector<Ciphertext> temp_product(YONGSOO_LIMIT);
  NTL_EXEC_RANGE(YONGSOO_LIMIT, first, last);
  for (int i = first; i < last; ++i) {
    evaluator.multiply_plain(B_cts_ntt[i], A_pts_ntt[i], temp_product[i]);
  }
  NTL_EXEC_RANGE_END

  Ciphertext temp_product_ntt = temp_product[0];
  for (size_t i = 1; i < YONGSOO_LIMIT; ++i) {
    evaluator.add_inplace(temp_product_ntt, temp_product[i]);
  }

  evaluator.transform_from_ntt(temp_product_ntt, C);
}

/*
 * Extending matrix multiplication of CCS18 to large blocks
 *
 * @param A:             vector of ciphertext A representing matrix of size
 * k1*k2. Each ciphertext has YONGSOO_LIMIT number of plaintext values. Matrix A
 * is values of \tilde{A} i.e., the output of the first transformation
 * completed. All such matrices are in row major blocks.
 * @param B:             vector of ciphertext B representing matrix of size
 * k2*k3
 * @param k1, k2, k3:    Matrix sizes
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix size of C is k1*k3
 */
void matrix_mult_general(const vector<Ciphertext> &A,
                         const vector<Ciphertext> &B, vector<Ciphertext> &C,
                         const size_t k1, const size_t k2, const size_t k3,
                         Evaluator &evaluator, BatchEncoder &batch_encoder,
                         GaloisKeys &gal_keys, RelinKeys &relin_keys) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

  NTL_EXEC_RANGE(k1, first, last);
  for (long rows = first; rows < last; ++rows) {
    Ciphertext temp;
    for (size_t columns = 0; columns < k3; ++columns) {
      for (size_t common_dim = 0; common_dim < k2; ++common_dim) {
        matrix_mult_unit(A[rows * k2 + common_dim],
                         B[common_dim * k3 + columns], temp, evaluator,
                         batch_encoder, gal_keys, relin_keys);
        evaluator.add_inplace(C[rows * k3 + columns], temp);
      }
    }
  }
  NTL_EXEC_RANGE_END;
}

/*
 * Extending matrix multiplication of CCS18 to large blocks
 *
 * @param A:             vector of plaintext A representing matrix of size
 * k1*k2. Each plaintext has YONGSOO_LIMIT number of data values. Matrix A is
 * values of \tilde{A} i.e., the output of the first transformation completed.
 * All such matrices are in row major blocks.
 * @param B:             vector of ciphertext B representing matrix of size
 * k2*k3
 * @param k1, k2, k3:    Matrix sizes
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix size of C is k1*k3
 */
void matrix_mult_general(const vector<Plaintext> &A,
                         const vector<Ciphertext> &B, vector<Ciphertext> &C,
                         const size_t k1, const size_t k2, const size_t k3,
                         Evaluator &evaluator, BatchEncoder &batch_encoder,
                         GaloisKeys &gal_keys) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

  NTL_EXEC_RANGE(k1, first, last);
  for (long rows = first; rows < last; ++rows) {
    Ciphertext temp;
    for (size_t columns = 0; columns < k3; ++columns) {
      for (size_t common_dim = 0; common_dim < k2; ++common_dim) {
        matrix_mult_unit(A[rows * k2 + common_dim],
                         B[common_dim * k3 + columns], temp, evaluator,
                         batch_encoder, gal_keys);
        evaluator.add_inplace(C[rows * k3 + columns], temp);
      }
    }
  }
  NTL_EXEC_RANGE_END
}

/*
 * Extending CCS18 matrix multiplication to large blocks in NTT space
 *
 * @param A:             vector of plaintext A representing matrix of size
 * k1*k2. Each plaintext has YONGSOO_LIMIT number of data values. Matrix A is
 * values of \tilde{A} i.e., the output of the first transformation completed.
 * All such matrices are in row major blocks.
 * @param B:             vector of ciphertext B representing matrix of size
 * k2*k3
 * @param k1, k2, k3:    Matrix sizes
 * @param evaluator:     HE Evaluator
 * @param batch_encoder: HE batch_encoder
 *
 * @return C:            vector of ciphertext C=AxB representing the product
 * matrix size of C is k1*k3
 */
void matrix_mult_general_ntt(const vector<Plaintext> &A,
                             const vector<Ciphertext> &B, vector<Ciphertext> &C,
                             const size_t k1, const size_t k2, const size_t k3,
                             Evaluator &evaluator, BatchEncoder &batch_encoder,
                             GaloisKeys &gal_keys) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

  NTL_EXEC_RANGE(k1, first, last);
  for (long rows = first; rows < last; ++rows) {
    Ciphertext temp;
    for (size_t columns = 0; columns < k3; ++columns) {
      for (size_t common_dim = 0; common_dim < k2; ++common_dim) {
        matrix_mult_ntt_unit(A[rows * k2 + common_dim],
                             B[common_dim * k3 + columns], temp, evaluator,
                             batch_encoder, gal_keys);
        evaluator.add_inplace(C[rows * k3 + columns], temp);
      }
    }
  }
  NTL_EXEC_RANGE_END
}

/*
 * Take inputs as all the plaintext/ciphertext components from the initial
 * ciphertexts of A (whole matrix) and B (one column) perform blockwise matrix
 * multiplication as A * B_col
 */
void matrix_mult_unit_colwise(const vector<vector<Plaintext>> &A_pts,
                              const vector<vector<Ciphertext>> &B_col_cts,
                              vector<Ciphertext> &C, Evaluator &evaluator) {
  size_t k1 = C.size();
  size_t k2 = A_pts.size() / k1;
  assert(k2 == B_col_cts.size());

  NTL_EXEC_RANGE(k1, first, last);
  for (int i = first; i < last; ++i) {
    Ciphertext temp_product;
    Ciphertext res;
    res.reserve(3);
    temp_product.reserve(3);

    // A[i * k2][-] * B[0][-]
    evaluator.multiply_plain(B_col_cts[0][0], A_pts[i * k2][0], res);
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
      evaluator.multiply_plain(B_col_cts[0][k], A_pts[i * k2][k], temp_product);
      evaluator.add_inplace(res, temp_product);
    }

    // A[i * k2 + j][-] * B[j][-]
    for (size_t j = 1; j < k2; ++j) {
      for (size_t k = 0; k < YONGSOO_LIMIT; ++k) {
        evaluator.multiply_plain(B_col_cts[j][k], A_pts[i * k2 + j][k],
                                 temp_product);
        evaluator.add_inplace(res, temp_product);
      }
    }

    C[i] = res;
  }
  NTL_EXEC_RANGE_END
}

/*
 * Take inputs as all the ciphertext components from the initial ciphertexts of
 * A (one row) and B (whole matrix) perform blockwise matrix multiplication as
 * A_row * B
 */
void matrix_mult_unit_rowwise(const vector<vector<Ciphertext>> &A_row_cts,
                              const vector<vector<Ciphertext>> &B_cts,
                              vector<Ciphertext> &C, Evaluator &evaluator,
                              RelinKeys &relin_keys) {
  size_t k3 = C.size();
  size_t k2 = B_cts.size() / k3;
  assert(k2 == A_row_cts.size());

  NTL_EXEC_RANGE(k3, first, last);
  for (int j = first; j < last; ++j) {
    Ciphertext temp_product;
    Ciphertext res;
    res.reserve(3);
    temp_product.reserve(3);

    evaluator.multiply(A_row_cts[0][0], B_cts[j][0], res);
    for (size_t k = 1; k < YONGSOO_LIMIT; ++k) {
      evaluator.multiply(A_row_cts[0][k], B_cts[j][k], temp_product);
      evaluator.add_inplace(res, temp_product);
    }

    for (size_t i = 1; i < k2; ++i) {
      for (size_t k = 0; k < YONGSOO_LIMIT; ++k) {
        evaluator.multiply(A_row_cts[i][k], B_cts[i * k3 + j][k], temp_product);
        evaluator.add_inplace(res, temp_product);
      }
    }
    evaluator.relinearize_inplace(res, relin_keys);
    C[j] = res;
  }
  NTL_EXEC_RANGE_END
}

void matrix_mult_general_ntt_fast(
    const vector<Plaintext> &A, const vector<Ciphertext> &B,
    vector<Ciphertext> &C, const size_t k1, const size_t k2, const size_t k3,
    Evaluator &evaluator, BatchEncoder &batch_encoder, GaloisKeys &gal_keys) {
  assert(A.size() == k1 * k2 && "Matrix A of incorrect dimensions");
  assert(B.size() == k2 * k3 && "Matrix B of incorrect dimensions");
  assert(C.size() == k1 * k3 && "Matrix C of incorrect dimensions");

  // generate all the ciphertext components of Bmat
  vector<vector<Ciphertext>> B_cts(k2 * k3, vector<Ciphertext>(YONGSOO_LIMIT));
  generate_all_ntt_components_Bcts(B_cts, B, evaluator, gal_keys);

  NTL_EXEC_RANGE(k1, first, last);
  for (int i = first; i < last; ++i) {
    // choose each column and store it as a vector of ciphertexts
    vector<Plaintext> A_row;
    for (size_t j = 0; j < k2; ++j) {
      A_row.push_back(A[i * k2 + j]);
    }
    vector<vector<Plaintext>> A_row_pts(k2, vector<Plaintext>(YONGSOO_LIMIT));
    generate_all_ntt_components_Apts(A_row_pts, A_row, B_cts[0][0], evaluator,
                                     batch_encoder);

    // perform matrix computation at one time
    vector<Ciphertext> res(k3);
    matrix_mult_unit_rowwise(A_row_pts, B_cts, res, evaluator);

    // store the computatio results in
    for (size_t j = 0; j < k3; ++j) {
      evaluator.transform_from_ntt(res[j], C[k3 * i + j]);
    }
  }
  NTL_EXEC_RANGE_END
}

void mult_wrapper(const size_t k1, const size_t k2, const size_t k3,
                  const size_t selection, const size_t logp) {
  size_t poly_modulus_degree = POLY_MODULUS_DEGREE;
  SUB_OPTIMAL = (POLY_MODULUS_DEGREE / (2 * YONGSOO_LIMIT * YONGSOO_LIMIT));
  assert(SUB_OPTIMAL == 1 or
         SUB_OPTIMAL == 2 && "Some issue with sub_optimality factor");

  // miran
  SetNumThreads(NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  string option;

  if (selection == 1)
    option = "ct-ct";
  else if (selection == 2)
    option = "pt-ct (NTT)";
  else if (selection == 3)
    option = "pt-ct (Standard)";

  string banner = option + " BFV matrix multiplication with Degree: ";

#if (HOISTING)
  option += " w/hoisting";
#endif
#if (FASTBLOCK)
  option += " w/fast-block-comp";
#endif

  print_example_banner(banner + to_string(poly_modulus_degree));

  int logn = log2(poly_modulus_degree);
  std::shared_ptr<SmallPContext> smallpcontext =
      setup_params(logp, logn, selection);

  // EncryptionParameters parms(scheme_type::BFV);
  // parms.set_poly_modulus_degree(poly_modulus_degree);
  // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  // parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, logp));
  // shared_ptr<SEALContext> context = SEALContext::Create(parms);

#if (VERBOSE)
  print_parameters(context);
  cout << endl;
#endif

  //   KeyGenerator keygen(context);
  //   auto secret_key = keygen.secret_key();
  //   auto public_key = keygen.public_key();

  //   Encryptor encryptor(context, public_key);
  //   Decryptor decryptor(context, secret_key);
  //   Evaluator evaluator(context);
  //   BatchEncoder batch_encoder(context);
  //   size_t slot_count = batch_encoder.slot_count();
  //   RelinKeys relin_keys = keygen.relin_keys();
  //   print_parameters(context);
  //   cout << endl;
  //   GaloisKeys gal_keys;

  //   /*****************************************************************
  //    Generate all the required Galois keys if GALKEY = 1;
  //    otherwise, generate power-of-two Galois keys as SEAL-default
  //    ******************************************************************/

  // #if (GALKEY)
  //   if (selection == 1) {
  //     vector<int> steps(3 * (YONGSOO_LIMIT - 1), 0);
  //     for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
  //       steps[i] = (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT;
  //       steps[i + (YONGSOO_LIMIT - 1)] = SUB_OPTIMAL * (-YONGSOO_LIMIT + i +
  //       1); steps[i + (YONGSOO_LIMIT - 1) * 2] = SUB_OPTIMAL * (i + 1);
  //     }
  //     gal_keys = keygen.galois_keys(steps);
  //   } else if ((selection == 2) || (selection == 3)) {
  //     vector<int> steps(YONGSOO_LIMIT - 1, 0);
  //     for (size_t i = 0; i < YONGSOO_LIMIT - 1; ++i) {
  //       steps[i] = (i + 1) * SUB_OPTIMAL * YONGSOO_LIMIT;
  //     }
  //     gal_keys = keygen.galois_keys(steps);
  //   }
  // #else
  // #if (HOISTING)
  //   throw invalid_argument("hoisting should prepare all the rotation keys");
  // #endif
  //   gal_keys = keygen.galois_keys();
  // #endif

  /*****************************************************************
   Generate sample matrices and the corresponding \tilde{} matrices
   Ensure A and B are valid plaintext (negative numbers wrapped
   around 2^64 but we want them wrapped around plain_modulus.value())
   Repeat B matrix for free rotations.
   ******************************************************************/

  size_t slot_count = smallpcontext->batch_encoder->slot_count();
  size_t d = 64;
  vector<uint64_t> data_1(slot_count, 0), A_data(d * d, 0),
      data_2(slot_count, 0), B_data(d * d, 0);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, 7);

  for (size_t i = 0; i < d * d; ++i) {
    //        A_data[i] = i;
    //        B_data[i] = d*d - i;
    A_data[i] = distribution(generator);
    B_data[i] = distribution(generator);
  }

  generate_tilde<uint64_t>(A_data, B_data, data_1, data_2, d, slot_count);

  /*****************************************************************
   Generate the initial plaintext as well as the encryptions of
   matrices \tilde{A} and \tilde{B}. Also generate an encryption of 0's
   *****************************************************************/
  Plaintext plain_1;
  // Plaintext plain_1(parms.poly_modulus_degree(), 0);
  smallpcontext->batch_encoder->encode(data_1, plain_1);

  Plaintext plain_2;
  // Plaintext plain_2(parms.poly_modulus_degree(), 0);
  smallpcontext->batch_encoder->encode(data_2, plain_2);

  Ciphertext encrypted_1(smallpcontext->seal_stuff->context);
  smallpcontext->seal_stuff_client->encryptor->encrypt(plain_1, encrypted_1);

  Ciphertext encrypted_2(smallpcontext->seal_stuff->context);
  smallpcontext->seal_stuff_client->encryptor->encrypt(plain_2, encrypted_2);
  cout << "Noise budget in fresh encrypted matrix: "
       << smallpcontext->seal_stuff_client->decryptor->invariant_noise_budget(
              encrypted_1)
       << " bits" << endl;
  cout << "--------------------------------------------" << endl;

  /*****************************************************************
   Generate the large matrices by pasting copies of the example
   matrices. Ensure matrix C is initialized as zeros.
   Also measure the time for this.
   *****************************************************************/
#if (VERBOSE)
  if (k1 < NUM_THREADS)
    cout << "Parallelization is sub-optimal. Running multiplication with " << k1
         << " threads" << endl;
  else
    cout << "Running multiplication with " << NUM_THREADS << " threads" << endl;
#endif
  cout << endl;

  chrono::high_resolution_clock::time_point time_start, time_end;
  chrono::TIME_UNIT time_total(0);

  vector<Ciphertext> C(k1 * k3, Ciphertext(smallpcontext->seal_stuff->context));

  auto &evaluator = *smallpcontext->seal_stuff->evaluator;
  auto &batch_encoder = *smallpcontext->batch_encoder;
  auto gal_keys = smallpcontext->seal_stuff_client->gal_keys;
  auto relin_keys = smallpcontext->seal_stuff_client->relin_keys;

  std::shared_ptr<BigPContext> bigpcontext = nullptr;

  if (selection == 1) {
    // RelinKeys relin_keys = smkeygen.relin_keys();
    vector<Ciphertext> A(k1 * k2, encrypted_1), B(k2 * k3, encrypted_2);
    time_start = chrono::high_resolution_clock::now();

#if (FASTBLOCK)
    matrix_mult_general_fast_ctct(A, B, C, k1, k2, k3, smallpcontext,
                                  bigpcontext);
#else
    matrix_mult_general(A, B, C, k1, k2, k3, evaluator, batch_encoder, gal_keys,
                        relin_keys);
#endif
  } else if (selection == 2) {
    vector<Plaintext> A(k1 * k2, plain_1);
    vector<Ciphertext> B(k2 * k3, encrypted_2);

    time_start = chrono::high_resolution_clock::now();

#if (FASTBLOCK)
    matrix_mult_general_ntt_fast(A, B, C, k1, k2, k3, evaluator, batch_encoder,
                                 gal_keys);
#else
    matrix_mult_general_ntt(A, B, C, k1, k2, k3, evaluator, batch_encoder,
                            gal_keys);
#endif
  }
  if (selection == 3) {
    vector<Plaintext> A(k1 * k2, plain_1);
    vector<Ciphertext> B(k2 * k3, encrypted_2);
    time_start = chrono::high_resolution_clock::now();

#if (FASTBLOCK)
    matrix_mult_general_fast_ptct(A, B, C, k1, k2, k3, smallpcontext);
#else
    matrix_mult_general(A, B, C, k1, k2, k3, evaluator, batch_encoder,
                        gal_keys);
#endif
  }

  time_end = chrono::high_resolution_clock::now();
  time_total = chrono::duration_cast<chrono::TIME_UNIT>(time_end - time_start);
  cout << "--------------------------------------------" << endl;
  cout << "Total time (" + option + "): \t\t" << time_total.count()
       << PRINT_STRING(TIME_UNIT) << endl;
  cout << "--------------------------------------------" << endl;

  /*****************************************************************
   Decrypt the matrix for C = A x B, unbatch it and check the answer
   *****************************************************************/
  Plaintext plain_3(poly_modulus_degree, 0);
  smallpcontext->seal_stuff_client->decryptor->decrypt(C[0], plain_3);
  cout << "Noise budget in product matrix: "
       << smallpcontext->seal_stuff_client->decryptor->invariant_noise_budget(
              C[0])
       << " bits" << endl;
  cout << "--------------------------------------------" << endl;

  vector<uint64_t> data_3(slot_count);
  batch_encoder.decode(plain_3, data_3);

  // decode a matrix
  vector<uint64_t> expected(d * d, 0);
  decode_matrix<uint64_t>(expected, data_3, d);

  // plain matrix
  vector<uint64_t> C_data(d * d, 0);
  mul<uint64_t>(A_data, B_data, C_data);
  mul<uint64_t>(C_data, C_data, k2);

  cout << "ok? " << check_equal<uint64_t>(C_data, expected) << endl;

#if (VERBOSE)
  cout << "Plain Comp" << endl;
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < d; ++j) cout << "\t" << C_data[i * d + j];
    cout << endl;
  }
  cout << "--------------------------------------------" << endl;
  cout << "HE Comp = " << k2 << "* Plain" << endl;
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < d; ++j)
      cout << "\t" << data_3[SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)];
    cout << endl;
  }
  cout << "        ------------" << endl;
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < d; ++j)
      cout << "\t"
           << data_3[slot_count - SUB_OPTIMAL * YONGSOO_LIMIT * d +
                     SUB_OPTIMAL * (i * YONGSOO_LIMIT + j)];
    cout << endl;
  }
  cout << "--------------------------------------------" << endl;
#endif
}
