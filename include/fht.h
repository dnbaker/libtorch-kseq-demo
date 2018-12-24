#include "vec/vec.h"

/*
def fwht(a):
     """In-place Fast Walsh-Hadamard Transform of array a"""
     h = 1
     while h < len(a):
        for i in range(0, len(a), h * 2):
           for j in range(i, i + h):
              x = a[j]
              y = a[j+h]
              a[j] = x + y
              a[j+h] = x - y
        h *= 2
*/

template<typename T>
void dumb_fht(T *buf, unsigned log_n) {
    size_t n = 1ull << log_n;
    for(unsigned i = 0; i < log_n; ++i) {
        size_t s1 = 1ull << i;
        size_t s2 = s1 << 1;
        for (int j = 0; j < n; j += s2) {
            for (int k = 0; k < s1; ++k) {
                auto u = buf[j + k];
                auto v = buf[j + k + s1];
                buf[j + k] = u + v;
                buf[j + k + s1] = u - v;
            }
        }
    }
}

#define DO_THING(k, s1, offset) do {\
            u = buf[j + k + offset];\
            v = buf[j + k + s1 + offset];\
            buf[j + k + offset] = u + v;\
            buf[j + k + s1 + offset] = u - v; }\
         while(0)

template<typename T>
void dumb_fht_first16(T *buf, unsigned log_n) {
    size_t n = size_t(1) << log_n;
    size_t s1 = 1ull << 3;
    size_t s2 = s1 << 1;
    for (int j = 0; j < n; j += s2) {
        {
            T u, v;
            DO_THING(0, 1, 0);
            DO_THING(1, 1, 0);
            DO_THING(0, 1, 2);
            DO_THING(1, 1, 2);
            DO_THING(0, 1, 4);
            DO_THING(1, 1, 4);
            DO_THING(0, 1, 6);
            DO_THING(1, 1, 6);
            DO_THING(0, 1, 8);
            DO_THING(1, 1, 8);
            DO_THING(0, 1, 10);
            DO_THING(1, 1, 10);
            DO_THING(0, 1, 12);
            DO_THING(1, 1, 12);
            DO_THING(0, 1, 14);
            DO_THING(1, 1, 14);
            DO_THING(0, 2, 0);
            DO_THING(1, 2, 0);
            DO_THING(2, 2, 0);
            DO_THING(3, 2, 0);
            DO_THING(0, 2, 4);
            DO_THING(1, 2, 4);
            DO_THING(2, 2, 4);
            DO_THING(3, 2, 4);
            DO_THING(0, 2, 8);
            DO_THING(1, 2, 8);
            DO_THING(2, 2, 8);
            DO_THING(3, 2, 8);
            DO_THING(0, 2, 12);
            DO_THING(1, 2, 12);
            DO_THING(2, 2, 12);
            DO_THING(3, 2, 12);
            DO_THING(0, 4, 0);
            DO_THING(1, 4, 0);
            DO_THING(2, 4, 0);
            DO_THING(3, 4, 0);
            DO_THING(4, 4, 0);
            DO_THING(5, 4, 0);
            DO_THING(6, 4, 0);
            DO_THING(7, 4, 0);
            DO_THING(0, 4, 8);
            DO_THING(1, 4, 8);
            DO_THING(2, 4, 8);
            DO_THING(3, 4, 8);
            DO_THING(4, 4, 8);
            DO_THING(5, 4, 8);
            DO_THING(6, 4, 8);
            DO_THING(7, 4, 8);
        }
    }
}

template<typename T>
void dumb_fht(T *buf, unsigned log_n) {
    size_t n = 1ull << log_n;
    for(unsigned i = 0; i < log_n; ++i) {
        size_t s1 = 1ull << i;
        size_t s2 = s1 << 1;
        for (int j = 0; j < n; j += s2) {
            for (int k = 0; k < s1; ++k) {
                auto u = buf[j + k];
                auto v = buf[j + k + s1];
                buf[j + k] = u + v;
                buf[j + k + s1] = u - v;
            }
        }
    }
}

template<typename IType>
void ffht_baseline(IType *p, unsigned l2sz) {
    assert(l2sz < 48 || !std::fprintf(stderr, "can't do sizes > 2**48, wow that's a lot of stuff."));
    size_t h = 1;
    const auto nelem = size_t(1) << l2sz;
    for(size_t h = 1, nelem = size_t(1) << l2sz; h != nelem; h <<= 1) {
        size_t h2 = h * 2;
        for(size_t i = 0; i < nelem; i += h2) {
            for(size_t j = i; j != i + h; ++j) {
                auto x = p[j], y = p[j+h];
                p[j] = x + y, p[j + h] = x - y;
            }
        }
    }
}
