#ifndef UTILS_H
#define UTILS_H

// Simple implementations of quicksort and binary search, using templates to
// have data type felixibility.

// Implementations are in this header file because they use templates.

/*============================================================================*/

/* Quicksort */

template<typename TYPE> void quicksort (TYPE[], int);
template<typename TYPE> void do_quicksort (TYPE[], int, int);
template<typename TYPE> int partition(TYPE[], int, int);
template<typename TYPE> int choose_pivot (TYPE[], int, int);

/* Binary search implementation */

template<typename TYPE> int binary_search (TYPE[], int, TYPE);
template<typename TYPE> int do_search (TYPE[], int, int, TYPE);

/*============================================================================*/

/*=========================\
| Quicksort implementation |
\=========================*/

// Main wrapper for Quicksort
// Receives a list and list size, and sorts the list in-place
template<typename TYPE>
void quicksort (TYPE data[], int size) {
  do_quicksort(data, 0, size-1);
}

// Recursively sorts the subarray between indices left and right, inclusive.
template<typename TYPE>
void do_quicksort (TYPE data[], int left, int right) {
  int pivot_idx;
  if (left < right) {
    pivot_idx = partition(data, left, right);
    do_quicksort(data, left, pivot_idx-1);
    do_quicksort(data, pivot_idx+1, right);
  }
}

// Partitions the subarray between indices left and right, inclusive.
// Returns the final index of the pivot.
template<typename TYPE>
int partition(TYPE data[], int left, int right) {
  int pivot_idx, next_idx, i;
  TYPE pivot_val, temp_val;
  pivot_idx = choose_pivot(data, left, right);
  pivot_val = data[pivot_idx];
  data[pivot_idx] = data[right];
  data[right] = pivot_val;
  next_idx = left;
  for (i=left; i<=right-1; i++) {
    if (data[i] <= pivot_val) {
      temp_val = data[next_idx];
      data[next_idx] = data[i];
      data[i] = temp_val;
      next_idx++;
    }
  }
  temp_val = data[next_idx];
  data[next_idx] = data[right];
  data[right] = temp_val;
  return next_idx;
}

// Chooses a pivot in the subarray from left to right, inclusive.
// Will choose the 'median of three' from the leftmost, central, and
// rightmost elements. Will also *sort* these three elements before returning.
// Returns the index of the chosen pivot.
template<typename TYPE>
int choose_pivot (TYPE data[], int left, int right) {
  int center;
  TYPE temp_val;
  center = (left+right)/2;
  if (data[left] > data[center]) {
    temp_val = data[center];
    data[center] = data[left];
    data[left] = temp_val;
  }
  if (data[left] > data[right]) {
    temp_val = data[right];
    data[right] = data[left];
    data[left] = temp_val;
  }
  if (data[center] > data[right]) {
    temp_val = data[right];
    data[right] = data[center];
    data[center] = temp_val;
  }
  return center;
}

/*============================================================================*/

/*=============================\
| Binary search implementation |
\=============================*/

// Main wrapper for binary search
// Receives an array, its size and an item, and searches for the item
// in the array. Returns the index of an occurrence of the item if it
// is found, -1 otherwise. Naturally, assumes the provided array is sorted.
// If multiple occurrences exist, the first found will be returned.
template<typename TYPE>
int binary_search (TYPE data[], int size, TYPE item) {
  return do_search(data, 0, size-1, item);
}

// Recursive binary search
// Searches for item between indices left and right, inclusive.
// Returns index of item if it is found, -1 otherwise.
template<typename TYPE>
int do_search (TYPE data[], int left, int right, TYPE item) {
  int center;
  if (left == right) {
    if (data[left] == item) return left;
    else return -1;
  }
  if (item < data[left] || item > data[right]) return -1;
  center = (left+right)/2;
  if (data[center] == item) return center;
  else if (data[center] > item) return do_search(data, left, center-1, item);
  else return do_search(data, center+1, right, item);
}

/*============================================================================*/

#endif // UTILS_H
