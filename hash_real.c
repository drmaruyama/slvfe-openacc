#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>

uint64_t hash_double(const double *v_, const size_t n)
{
  uint64_t x = 0, buf;
  size_t i;
  for(i = 0; i < n; ++i){
    memcpy(&buf, &v_[i], sizeof(double));
    x = ((x << 7) | (x >> 57)) ^ buf;
  }
  return x;
}

uint64_t hash_float(const float *v_, const size_t n)
{
  uint64_t x = 0;
  uint32_t buf;
  size_t i;
  for(i = 0; i < n; ++i){
    memcpy(&buf, &v_[i], sizeof(float));
    x = ((x << 7) | (x >> 57)) ^ (uint64_t) buf;
  }
  return x;
}
