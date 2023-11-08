#include "quicksort.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "get_time.h"

using Type = long long;

inline uint64_t hash64(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

int main(int argc, char* argv[]) {
  size_t n = 1e8;
  int num_rounds = 3;
  int random_seed = 1;
  int test = 1;
  if (argc >= 2) {
    test = atoi(argv[1]);
  }
  if (argc >= 3) {
    n = atoll(argv[2]);
  }
  if (argc >= 4) {
    num_rounds = atoi(argv[3]);
  }
  if (argc >= 5) {
    random_seed = atoi(argv[4]);
  }

  Type* A = (Type*)malloc(n * sizeof(Type));
  Type* B = (Type*)malloc(n * sizeof(Type));

  double total_time = 0;
  for (int i = 0; i <= num_rounds; i++) {
    if(test == 2) {
      // Generate random arrays: test2
      parallel_for(0, n,
                 [&](size_t j) { A[j] = B[j] = hash64(j * random_seed) % (long unsigned int)1e5; });
    } else if(test == 4) {
      // Generate random arrays: test4
      parallel_for(0, n,
                   [&](size_t j) { A[j] = B[j] = hash64(j * random_seed) % (long unsigned int)10; });
    } else {
      // Generate random arrays: test1
      parallel_for(0, n,
                  [&](size_t j) { A[j] = B[j] = hash64(j * random_seed); });
    }

    // for(size_t j = 0; j < n; j++ ) {
    //   std::cout << A[j] << std::endl;
    // }
    parlay::timer t;
    quicksort(A, n);
    t.stop();

    std::sort(B, B + n);

    // std::cout << "Manual sort, Auto Sort" << std::endl;
    // for(size_t j = 0; j < n; j++ ) {
    //   std::cout << A[j] << " " << B[j] << std::endl;
    // }

    parallel_for(0, n, [&](size_t j) {
      if (A[j] != B[j]) {
        std::cout << "The output is not sorted\n";
        exit(0);
      }
    });

    if (i == 0) {
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time += t.total_time();
    }
  }
  std::cout << "Average running time: " << total_time / num_rounds << std::endl;

  free(A);
  return 0;
}
