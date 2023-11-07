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

template <class T>
void sequential_quicksort(T *A, size_t n) {
  if(n <= 1) {
    return;
  }
  T temp;
  size_t random_index = rand() % n;
  T pivot = A[random_index];
  size_t pivot_index = -1;
  for (size_t i = 0; i < n; i++) {
    if (A[i] <= pivot)
      pivot_index++;
  }
  // Swap the pivot to pivot index
  A[random_index] = A[pivot_index];
  A[pivot_index] = pivot;

  // Sorting left and right part of the pivot element
  size_t i = 0, j = n-1;
  while(i < pivot_index && j > pivot_index) {
    while(A[i] <= pivot) {
      i++;
    }
    while(A[j] > pivot) {
      j--;
    }
    if(i < pivot_index && j > pivot_index) {
      // Swap the numbers from left part to right part of the pivot table
      temp = A[i];
      A[i] = A[j];
      A[j] = temp;
      i++;
      j--;
    }
  }
  // std::cout << "indices " << 0 << " " << pivot_index << " " << n-1 << std::endl;
  // for(size_t i = 0; i<n; i++) {
  //   std::cout << A[i] << " ";
  // }
  // std::cout << std::endl;
  sequential_quicksort(A, pivot_index);
  sequential_quicksort(A + pivot_index + 1, n - pivot_index - 1);
}

// 2.68952
template <typename T>
T scan3(T *A, size_t n) {
  size_t k = sqrt(n);
  size_t chunk_size = (size_t) n/k + 1;
  T *chunk_sum = (T *)malloc(k * sizeof(T));
  T total = 0;
  T temp, end_of_chunk;

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

  free(chunk_sum);

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

template <typename T>
T scan_up(T* A, T* LS, size_t n) {
  if(n <= 3e6) {
    T sum = 0;
    for(size_t i = 0; i < n; i++) {
      sum += A[i];
    }
    return sum;
  }
  size_t m = n/2;
  T l,r;
  auto f1 = [&]() { l = scan_up(A, LS, m); };
  auto f2 = [&]() { r = scan_up(A + m, LS + m, n - m); };
  par_do(f1, f2);

  LS[m-1] = l;
  return l+r;
}

template <typename T>
void scan_down(T* A, T* LS, size_t n, T offset) {
  if(n <= 3e6) {
    T total = offset;
    for (size_t i = 0; i < n; i++) {
      T tmp = A[i];
      A[i] = total;
      total += tmp;
    }
    return;
  }
  size_t m = n/2;
  auto f1 = [&]() { scan_down(A, LS, m, offset); };
  auto f2 = [&]() { scan_down(A + m, LS + m, n - m, offset + LS[m-1]); };
  par_do(f1, f2);
}

// 2.734
template <typename T>
T scan2(T *A, size_t n, T* LS) {
  T offset = 0;
  T sum = scan_up(A, LS, n);
  scan_down(A, LS, n, offset);
  return sum;
}

// 50: 
// 100: 1.65716
// 200: 
template <class T>
size_t find_median_index(T* A, size_t n) {
  size_t k = 200;
  size_t random_ind[k];
  for(size_t i = 0; i < k; i++) {
    random_ind[i] = hash164(i) % n;
  }
  // std::cout << "Before sorting: " << std::endl;
  // for(size_t i = 0; i < k; i++) {
  //   std::cout << A[random_ind[i]] << " ";
  // }
  // std::cout << std::endl;
  std::sort(random_ind, random_ind + k, [&](size_t a, size_t b) {
    return A[a] < A[b];
  });
  // std::cout << "After sorting: " << std::endl;
  // for(size_t i = 0; i < k; i++) {
  //   std::cout << A[random_ind[i]] << " ";
  // }
  // std::cout << std::endl;
  return random_ind[k/2];
}

template <class T>
size_t parallel_partition(T *A, size_t n, T* B, size_t* LSl, size_t* LSr, size_t* left_prefix_sum, size_t* right_prefix_sum) {
  size_t random_index = find_median_index(A, n);
  T pivot = A[random_index];
  size_t pivot_index;

  parallel_for(0, n, [&](size_t i) {
    left_prefix_sum[i] = (A[i] <= pivot && i != random_index) ? 1 : 0;
    right_prefix_sum[i] = A[i] > pivot ? 1 : 0;
    B[i] = A[i];
  });

  left_prefix_sum[n] = 0;
  right_prefix_sum[n] = 0;

  auto f1 = [&]() { pivot_index = scan2(left_prefix_sum, n + 1, LSl); };
  auto f2 = [&]() { scan2(right_prefix_sum, n + 1, LSr);};
  par_do(f1, f2);  

  // std::cout << "Random Index: " << random_index << std::endl;
  // std::cout << "Left Prefix, Right Prefix" << std::endl;
  // for(size_t j = 0; j < n; j++ ) {
  //   std::cout << left_prefix_sum[j] << " " << right_prefix_sum[j] << std::endl;
  // }
  parallel_for(0, n, [&](size_t i) {
    if(left_prefix_sum[i + 1] != left_prefix_sum[i]) {
      A[left_prefix_sum[i]] = B[i];
    }
    if(right_prefix_sum[i + 1] != right_prefix_sum[i]) {
      A[pivot_index + right_prefix_sum[i] + 1] = B[i];
    }
  });
  A[pivot_index] = pivot;
  return pivot_index;
}

// 3e5: 1.62398
// 35e4: 1.63964
// 4e5: 1.61909
// 5e5: 1.66368
template <class T>
void quicksort_rec(T *A, size_t n, T *B, size_t* LSl, size_t* LSr, size_t* left_prefix_sum, size_t* right_prefix_sum) {
  // std::sort(A, A + n);
  if(n <= 1) {
    return;
  }
  if(n < 400000) {
    // sequential_quicksort(A, n);
    std::sort(A, A + n);
    return;
  }
  size_t pivot_index = parallel_partition(A, n, B, LSl, LSr, left_prefix_sum, right_prefix_sum);
  // std::cout << "Pivot Index: " << pivot_index << std::endl;
  // std::cout << "After Partitioning:" << std::endl;
  // for(size_t j = 0; j < n; j++ ) {
  //   std::cout << A[j] << std::endl;
  // }
  auto f1 = [&]() { quicksort_rec(A, pivot_index, B, LSl, LSr, left_prefix_sum, right_prefix_sum); };
  auto f2 = [&]() { quicksort_rec(A + pivot_index + 1, n - pivot_index - 1, B + pivot_index + 1, LSl + pivot_index + 1, LSr + pivot_index + 1, left_prefix_sum + pivot_index + 1, right_prefix_sum + pivot_index + 1); };
  par_do(f1, f2);
}

template <class T>
void quicksort(T *A, size_t n) {
  T* B = (T*)malloc(n * sizeof(T));
  size_t *LSl = (size_t *)malloc((n-1) * sizeof(size_t));
  size_t *LSr = (size_t *)malloc((n-1) * sizeof(size_t));
  size_t* left_prefix_sum = (size_t*)malloc((n+1) * sizeof(size_t));
  size_t* right_prefix_sum = (size_t*)malloc((n+1) * sizeof(size_t));

  quicksort_rec(A, n, B, LSl, LSr, left_prefix_sum, right_prefix_sum);

  free(B);
  free(LSl);
  free(LSr);
  free(left_prefix_sum);
  free(right_prefix_sum);
}
