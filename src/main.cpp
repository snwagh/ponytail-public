/*******************************************************
 * 2019 Sameer Wagh <snwagh@gmail.com>
 *
 * Matrix multiplication using SEAL and BFV implemented
 * in this file.
 *
 * This file is part of ponytail project.
 *******************************************************/

#include "matmul.h"
#include "seal_augment.h"

using namespace std;
using namespace seal;

extern size_t NUM_THREADS;
void print_usage(const char *bin);

int main(int argc, char **argv) {
#ifdef SEAL_VERSION
  if (VERBOSE) cout << "Microsoft SEAL version: " << SEAL_VERSION << endl;
#endif
  //test_multiply(); 
  // test_encode_decode();
  // test_ntt(128, 15);
  //test_encrypt_decrypt();
  // return 0;
  //test_linear_monomials(32, 16);
  //return 0;
  if (argc < 7) print_usage(argv[0]);

  NUM_THREADS = atoi(argv[1]);
  size_t selection = atoi(argv[2]);
  size_t k1 = atoi(argv[3]);
  size_t k2 = atoi(argv[4]);
  size_t k3 = atoi(argv[5]);
  size_t logp = atoi(argv[6]);

  if (logp <= 60) {
    mult_wrapper(k1, k2, k3, selection, logp);
  } else {
    mult_wrapper_bigp(k1, k2, k3, selection, logp);
  }
  return 0;
}

void print_usage(const char *bin) {
  cout << "Usage: " << bin << " NUM_THREADS MULT_TYPE k1 k2 k3 k4" << endl;
  cout << endl;
  cout << "Required Arguments:\n";
  cout << "NUM_THREADS        \tNumber of threads to be used\n";
  cout << "MULT_TYPE          \t1 for Ct-Ct, 2 for Pt-Ct (NTT), 3 for Pt-Ct "
          "(Standard)\n";
  cout << "k1                 \tRows \n";
  cout << "k2                 \tCommon_dimension\n";
  cout << "k3                 \tColumns\n";
  cout << "k4                 \tlog(p)\n";

  cout << endl;
  cout << "Report bugs to swagh@princeton.edu" << endl;
  exit(-1);
}
