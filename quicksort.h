#include <algorithm>
#include <math.h>

#include "parallel.h"

using namespace parlay;

inline uint64_t hash164(uint64_t u) {
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

size_t scan(size_t *A, size_t n, size_t *chunk_sum) {
  size_t k = sqrt(n);
  size_t chunk_size = n/k + 1;

  size_t total = 0;
  size_t temp, end_of_chunk;

  parallel_for(0, k, [&](size_t i) {
    chunk_sum[i] = 0;
    for(size_t j = 0; j < chunk_size - 1 && (i * chunk_size + j) < n - 1; j++) {
      chunk_sum[i] += A[i * chunk_size + j];
    }
  });

  for(size_t i = 0; i < k; i++) {
    if(i + 1 == k) {
      end_of_chunk = n - 1;
    } else {
      end_of_chunk = (i + 1) * chunk_size - 1;
    }
    total += chunk_sum[i];
    temp = A[end_of_chunk];
    A[end_of_chunk] = total;
    total += temp;
  }

  parallel_for(0, k, [&](size_t i) {
    if(i + 1 == k) {
      end_of_chunk = n - 1;
    } else {
      end_of_chunk = (i + 1) * chunk_size - 1;
    }
    for(size_t j = end_of_chunk; j > i * chunk_size; j--) {
      A[j - 1] = A[j] - A[j - 1];
    }
  });

  return total;
}

template <class T>
size_t find_median_index(T* A, size_t n) {
  size_t k = 20;
  size_t random_ind[k];
  for(size_t i = 0; i < k; i++) {
    random_ind[i] = hash164(i) % n;
  }

  std::sort(random_ind, random_ind + k, [&](size_t a, size_t b) {
    return A[a] < A[b];
  });

  return random_ind[k/2 - 1];
}

template <class T>
void parallel_partition(
  T *A,
  size_t n,
  size_t *pivot_inds,
  T* B,
  size_t* left_chunk_sum,
  size_t* right_chunk_sum,
  size_t* left_prefix_sum,
  size_t* right_prefix_sum) {
  size_t random_index = find_median_index(A, n);
  T pivot = A[random_index];

  parallel_for(0, n, [&](size_t i) {
    left_prefix_sum[i] = (A[i] < pivot) ? 1 : 0;
    right_prefix_sum[i] = A[i] > pivot ? 1 : 0;
    B[i] = A[i];
  });

  left_prefix_sum[n] = 0;
  right_prefix_sum[n] = 0;
  auto f1 = [&]() { pivot_inds[0] = scan(left_prefix_sum, n + 1, left_chunk_sum); };
  auto f2 = [&]() { pivot_inds[1] = scan(right_prefix_sum, n + 1, right_chunk_sum);};
  par_do(f1, f2);  

  pivot_inds[1] = n - pivot_inds[1];

  parallel_for(0, n, [&](size_t i) {
    if(left_prefix_sum[i + 1] != left_prefix_sum[i]) {
      A[left_prefix_sum[i]] = B[i];
    } else if(right_prefix_sum[i + 1] != right_prefix_sum[i]) {
      A[pivot_inds[1] + right_prefix_sum[i]] = B[i];
    }
  });  
  parallel_for(pivot_inds[0], pivot_inds[1], [&](size_t i) {
    A[i] = pivot;
  });
}

template <class T>
void quicksort_rec(
  T *A,
  size_t n,
  T *B,
  size_t* left_chunk_sum,
  size_t* right_chunk_sum,
  size_t* left_prefix_sum,
  size_t* right_prefix_sum) {
  if(n < 8750000) {
    std::sort(A, A + n);
    return;
  }
  size_t pivot_inds[2];
  parallel_partition(
    A,
    n,
    pivot_inds,
    B,
    left_chunk_sum,
    right_chunk_sum,
    left_prefix_sum,
    right_prefix_sum);
  auto f1 = [&]() { quicksort_rec(
    A,
    pivot_inds[0],
    B,
    left_chunk_sum,
    right_chunk_sum,
    left_prefix_sum,
    right_prefix_sum); };
  auto f2 = [&]() { quicksort_rec(
    A + pivot_inds[1],
    n - pivot_inds[1],
    B + pivot_inds[1],
    left_chunk_sum + pivot_inds[1],
    right_chunk_sum + pivot_inds[1],
    left_prefix_sum + pivot_inds[1],
    right_prefix_sum + pivot_inds[1]); };
  par_do(f1, f2);
}

template <class T>
void quicksort(T *A, size_t n) {
  T* B = (T*)malloc(n * sizeof(T));
  size_t *left_chunk_sum = (size_t *)malloc((n - 1) * sizeof(size_t));
  size_t *right_chunk_sum = (size_t *)malloc((n - 1) * sizeof(size_t));
  size_t* left_prefix_sum = (size_t*)malloc((n+1) * sizeof(size_t));
  size_t* right_prefix_sum = (size_t*)malloc((n+1) * sizeof(size_t));

  quicksort_rec(
    A,
    n,
    B,
    left_chunk_sum,
    right_chunk_sum,
    left_prefix_sum,
    right_prefix_sum);

  free(B);
  free(left_chunk_sum);
  free(right_chunk_sum);
  free(left_prefix_sum);
  free(right_prefix_sum);
}
