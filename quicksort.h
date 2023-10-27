#include <algorithm>

#include "parallel.h"

using namespace parlay;

template <class T>
void quicksort(T *A, size_t n) {
  std::sort(A, A + n);
}
