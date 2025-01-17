#ifndef ZFP_CACHE2_HPP
#define ZFP_CACHE2_HPP

#include "zfp/internal/array/cache.hpp"
#include <ctime>
#include <iostream>
#include <ratio>
#include <chrono>
#if defined(__x86_64__) || defined(_M_X64)
#include <x86intrin.h>
#endif


double total_time_compress = 0;
double total_time_decompress = 0;
size_t total_decompression = 0, total_compression = 0;
uint64_t cycles_compression = 0, cycles_decompression = 0;

#if defined(__x86_64__) || defined(_M_X64)
uint64_t inline rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#endif

#if defined(__GNUC__) && defined(__riscv)  && !defined(HAVE_TICK_COUNTER)
static __inline__ uint64_t rdtsc(void)
{
        uint64_t cycles;
        __asm__ volatile("rdcycle %0" : "=r"(cycles));
        return cycles;
}
#endif

namespace zfp {
namespace internal {

template <typename Scalar, class Store>
class BlockCache2 {
public:
  // constructor of cache of given size
  BlockCache2(Store& store, size_t bytes = 0) :
    cache(lines(bytes, store.blocks())),
    store(store)
  {}

  // byte size of cache data structure components indicated by mask
  size_t size_bytes(uint mask = ZFP_DATA_ALL) const
  {
    size_t size = 0;
    size += cache.size_bytes(mask);
    if (mask & ZFP_DATA_META)
      size += sizeof(*this);
    return size;
  }

  // cache size in number of bytes (cache line payload data only)
  size_t size() const { return cache.size() * sizeof(CacheLine); }

  // set minimum cache size in bytes (inferred from blocks if zero)
  void resize(size_t bytes)
  {
    flush();
    cache.resize(lines(bytes, store.blocks()));
  }

  // empty cache without compressing modified cached blocks
  void clear() const { cache.clear(); }

  // flush cache by compressing all modified cached blocks
  void flush() const
  {
    for (typename zfp::internal::Cache<CacheLine>::const_iterator p = cache.first(); p; p++) {
      if (p->tag.dirty()) {
        size_t block_index = p->tag.index() - 1;
        store.encode(block_index, p->line->data());
      }
      cache.flush(p->line);
    }
  }

  // perform a deep copy
  void deep_copy(const BlockCache2& c) { cache = c.cache; }

  // inspector
  Scalar get(size_t i, size_t j) const
  {
    const CacheLine* p = line(i, j, false);
    return (*p)(i, j);
  }

  // mutator
  void set(size_t i, size_t j, Scalar val)
  {
    CacheLine* p = line(i, j, true);
    (*p)(i, j) = val;
  }

  // reference to cached element
  Scalar& ref(size_t i, size_t j)
  {
    CacheLine* p = line(i, j, true);
    return (*p)(i, j);
  }

  // read-no-allocate: copy block from cache on hit, else from store without caching
  void get_block(size_t block_index, Scalar* p, ptrdiff_t sx, ptrdiff_t sy) const
  {
    const CacheLine* line = cache.lookup((uint)block_index + 1, false);
    if (line)
      line->get(p, sx, sy, store.block_shape(block_index));
    else {
      uint64_t c1 = rdtsc();
      store.decode(block_index, p, sx, sy);
      uint64_t c2 = rdtsc();
      fprintf(stderr, "get %llu\n", c2-c1);
    }
  }

  // write-no-allocate: copy block to cache on hit, else to store without caching
  void put_block(size_t block_index, const Scalar* p, ptrdiff_t sx, ptrdiff_t sy)
  {
    CacheLine* line = cache.lookup((uint)block_index + 1, true);
    if (line)
      line->put(p, sx, sy, store.block_shape(block_index));
    else {
      uint64_t c1 = rdtsc();
      store.encode(block_index, p, sx, sy);
      uint64_t c2 = rdtsc();
      fprintf(stderr, "put %llu\n", c2-c1);
    }
  }

protected:
  // cache line representing one block of decompressed values
  class CacheLine {
  public:
    // accessors
    Scalar operator()(size_t i, size_t j) const { return a[index(i, j)]; }
    Scalar& operator()(size_t i, size_t j) { return a[index(i, j)]; }

    // pointer to decompressed block data
    const Scalar* data() const { return a; }
    Scalar* data() { return a; }

    // copy whole block from cache line
    void get(Scalar* p, ptrdiff_t sx, ptrdiff_t sy) const
    {
      const Scalar* q = a;
      for (uint y = 0; y < 4; y++, p += sy - 4 * sx)
        for (uint x = 0; x < 4; x++, p += sx, q++)
          *p = *q;
    }

    // copy partial block from cache line
    void get(Scalar* p, ptrdiff_t sx, ptrdiff_t sy, uint shape) const
    {
      if (!shape)
        get(p, sx, sy);
      else {
        // determine block dimensions
        uint nx = 4 - (shape & 3u); shape >>= 2;
        uint ny = 4 - (shape & 3u); shape >>= 2;
        const Scalar* q = a;
        for (uint y = 0; y < ny; y++, p += sy - (ptrdiff_t)nx * sx, q += 4 - nx)
          for (uint x = 0; x < nx; x++, p += sx, q++)
            *p = *q;
      }
    }

    // copy whole block to cache line
    void put(const Scalar* p, ptrdiff_t sx, ptrdiff_t sy)
    {
      Scalar* q = a;
      for (uint y = 0; y < 4; y++, p += sy - 4 * sx)
        for (uint x = 0; x < 4; x++, p += sx, q++)
          *q = *p;
    }

    // copy partial block to cache line
    void put(const Scalar* p, ptrdiff_t sx, ptrdiff_t sy, uint shape)
    {
      if (!shape)
        put(p, sx, sy);
      else {
        // determine block dimensions
        uint nx = 4 - (shape & 3u); shape >>= 2;
        uint ny = 4 - (shape & 3u); shape >>= 2;
        Scalar* q = a;
        for (uint y = 0; y < ny; y++, p += sy - (ptrdiff_t)nx * sx, q += 4 - nx)
          for (uint x = 0; x < nx; x++, p += sx, q++)
            *q = *p;
      }
    }

  protected:
    static size_t index(size_t i, size_t j) { return (i & 3u) + 4 * (j & 3u); }
    Scalar a[4 * 4];
  };

  // return cache line for (i, j); may require write-back and fetch
  CacheLine* line(size_t i, size_t j, bool write) const
  {
    // using namespace std::chrono;
    CacheLine* p = 0;
    size_t block_index = store.block_index(i, j);
    typename zfp::internal::Cache<CacheLine>::Tag tag = cache.access(p, (uint)block_index + 1, write);
    size_t stored_block_index = tag.index() - 1;
    uint64_t c1 = 0, c2 = 0;
    if (stored_block_index != block_index) {
      // write back occupied cache line if it is dirty
      if (tag.dirty()) {
        // high_resolution_clock::time_point t1 = high_resolution_clock::now();
        c1 = rdtsc();
        store.encode(stored_block_index, p->data());
        // high_resolution_clock::time_point t2 = high_resolution_clock::now();
        c2 = rdtsc();
        // duration<double> time_span_compress = duration_cast<duration<double> >(t2 - t1);
        // total_time_compress += time_span_compress.count();
        total_compression += 1;
        cycles_compression += c2 - c1;
        // std::cout << "compression cycles "<< c2-c1 << std::endl;
      }
      // fetch cache line
      // high_resolution_clock::time_point t3 = high_resolution_clock::now();
      uint64_t c3 = rdtsc();
      store.decode(block_index, p->data());
      uint64_t c4 = rdtsc();
      // high_resolution_clock::time_point t4 = high_resolution_clock::now();
      // duration<double> time_span_decompress = duration_cast<duration<double> >(t4 - t3);
      // total_time_decompress += time_span_decompress.count();
      total_decompression += 1;
      cycles_decompression += c4 - c3;
      fprintf(stderr, "%llu,%llu\n", c2-c1 , c4-c3);
    }
    return p;
  }

  // default number of cache lines for array with given number of blocks
  static uint lines(size_t blocks)
  {
    // compute m = O(sqrt(n))
    size_t m;
    for (m = 1; m * m < blocks; m *= 2);
    return static_cast<uint>(m);
  }

  // number of cache lines corresponding to size (or suggested size if zero)
  static uint lines(size_t bytes, size_t blocks)
  {
    // ensure block index fits in tag
    if (blocks >> ((sizeof(uint) * CHAR_BIT) - 1))
      throw zfp::exception("zfp array too large for cache");
    uint n = bytes ? static_cast<uint>((bytes + sizeof(CacheLine) - 1) / sizeof(CacheLine)) : lines(blocks);
    return std::max(n, 1u);
  }

  mutable Cache<CacheLine> cache; // cache of decompressed blocks
  Store& store;                   // store backed by cache
};

} // internal
} // zfp

#endif
