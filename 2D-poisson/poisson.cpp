// Solve Poisson equation div(grad(u)) = f(x, y) = sqrt(x^2 + y^2)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iostream>
#include <sys/time.h>
// #include "numrep.hpp"

#define GHOST 8 // 8 layers to accommodate Blaz

// grid size (excluding ghost layers)
#ifndef SIZE
  #define SIZE 1024
#endif

// number of iterations
#ifndef ITERS
  #define ITERS 32768
#endif

//#define JACOBI 1 // use Jacobi update when defined

// available number representations; set via Makefile
//#define REAL NumRep<unsigned short, ExpGolombRice<2>, MapLinear> // posit-16
//#define REAL NumRep<unsigned int, ExpGolombRice<2>, MapLinear>   // posit-32
//#define REAL __fp16                                              // IEEE half precision
//#define REAL float                                               // IEEE single precision
//#define REAL double                                              // IEEE double precision
//#define ZFP_RATE 12                                              // zfp rate in bits/value when defined
//#define BLAZ 1                                                   // use Blaz matrix compression when defined

#ifdef REAL
  typedef REAL real;
#elif defined(ZFP_RATE)
  #include "zfp/array2.hpp"
#elif BLAZ
  extern "C" {
    #include "blaz.h"
  }
#else
  // default: double precision
  #define REAL double
  typedef REAL real;
#endif

// timer
inline double
now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

// return u(i, j)
inline double
get(const double* u, size_t n, size_t i, size_t j)
{
  return u[i + n * j];
}

// set u(i, j) = val
template <typename real>
inline void
set(real* u, size_t n, size_t i, size_t j, real val)
{
  u[i + n * j] = val;
}

// compress/decompress u
inline void
store(double* u, size_t n, double* v = 0)
{
  if (!v)
    v = u;

#ifdef REAL
  // convert u to then from specified reduced-precision representation 'real'
  for (size_t i = 0; i < n * n; i++)
    u[i] = double(real(v[i]));
#elif defined(ZFP_RATE)
  // use zfp fixed-rate arrays to compress then decompress u
  zfp::array2<double> a(n, n, ZFP_RATE, v);
//  fprintf(stderr, "zfp rate = %g\n", a.rate());
  a.flush_cache();
  a.get(u);
#elif BLAZ
  // construct uncompressed Blaz matrix
  Blaz_Matrix* mat = (Blaz_Matrix*)blaz_malloc(sizeof(Blaz_Matrix));
  mat->width = n;
  mat->height = n;
  mat->matrix = v;
  // compress matrix
  Blaz_Compressed_Matrix* cmat = blaz_compress(mat);
  free(mat);
  // decompress matrix
  mat = blaz_uncompress(cmat);
  free(cmat->block_first_elts);
  free(cmat->block_mean_slope);
  free(cmat->compressed_values);
  free(cmat);
  // copy decompressed matrix to u
  std::copy(mat->matrix, mat->matrix + n * n, u);
  // clean up
  free(mat->matrix);
  free(mat);
#else
  #error "no number type defined"
#endif
}

// output gradient norm or Laplacian in single precision
void
write(const double* u, size_t m, size_t g, size_t frame)
{
  const size_t n = g + m + g; // number of grid cells including ghost layers
  const double h = 2. / m; // grid spacing for domain [-1, 1]^3

  // compute field, v
  float* v = new float[m * m];
  double error = 0;
  for (size_t j = g; j < n - g; j++) {
    double y = h * (j - g - 0.5 * (m - 1));
    for (size_t i = g; i < n - g; i++) {
      double x = h * (i - g - 0.5 * (m - 1));
      double r = std::sqrt(x * x + y * y);
#if FUNCTION
      // compute (9 u)^(1/3) = r
      double val = std::cbrt(9 * get(u, n, i, j));
#elif GRADIENT
      // compute gradient norm v = sqrt(3 |grad(u)|) = r
      double ux = (get(u, n, i + 1, j) - get(u, n, i - 1, j)) / (2 * h);
      double uy = (get(u, n, i, j + 1) - get(u, n, i, j - 1)) / (2 * h);
      double val = std::sqrt(3 * std::sqrt(ux * ux + uy * uy));
#else
      // compute Laplacian v = div(grad(u)) = r and solution error
      double uxx = get(u, n, i - 1, j) - 2 * get(u, n, i, j) + get(u, n, i + 1, j);
      double uyy = get(u, n, i, j - 1) - 2 * get(u, n, i, j) + get(u, n, i, j + 1);
      double val = (uxx + uyy) / (h * h);
#endif
      set(v, m, i - g, j - g, float(val));
      error += (val - r) * (val - r);
    }
  }
  error = std::sqrt(error / (m * m));
#if GRADIENT
  fprintf(stderr, "gradient error = %g\n", error);
#else
  fprintf(stderr, "L2 Laplacian error = %g\n", error);
#endif

  // output computed field, v, without ghost layers
  char path[0x1000];
  snprintf(path, 20, "/tmp/frame%04u.bof", (unsigned int)frame);
  FILE* file = fopen(path, "wb");
  if (file) {
    fwrite(v, sizeof(*v), m * m, file);
    fclose(file);
  }

  delete[] v;
}

inline void
symmetrize(double* u, size_t m, int g, size_t i, size_t j)
{
  const size_t n = g + m + g; // number of grid cells including ghost layers
  double val = u[i + n * j];
  u[(n - 1 - i) + n * j] = val;
  u[i + n * (n - 1 - j)] = val;
  u[(n - 1 - i) + n * (n - 1 - j)] = val;
}

// allocate and initialize state
double*
init(size_t m, int g)
{
  const size_t n = g + m + g; // number of grid cells including ghost layers
  const double h = 2. / m; // grid spacing for domain [-1, 1]^3

  double* u = new double[n * n];

  for (int j = 0; j < n; j++) {
    double y = h * (j - g - 0.5 * (m - 1));
    for (int i = 0; i < n; i++) {
      double x = h * (i - g - 0.5 * (m - 1));
      double r = std::sqrt(x * x + y * y);
#if 0
      // initialize with constant everywhere
      double val = 1. / 9;
#elif 0
      // initialize with final solution
      double val = r * r * r / 9;
#elif 0
      // initialize u(x, y) = r^3 / 9 on boundary, zero in interior
      double val = (i == g || i == n - 1 - g || j == g || j == n - 1 - g) ? r * r * r / 9 : 0.0;
#elif 0
      // initialize with interpolation of boundary
      double rx = std::sqrt(x * x + 1);
      double ry = std::sqrt(1 + y * y);
      double ux = rx * rx * rx / 9;
      double uy = ry * ry * ry / 9;
      double val = ((1 - x * x) * (4 - x * x + 5 * y * y) * ux +
                    (1 - y * y) * (4 - y * y + 5 * x * x) * uy) /
                   (16 - 8 * (x * x + y * y) + (x * x - y * y) * (x * x - y * y));
#elif 0
      // blend true solution, f, with another function, g
      double f = r * r * r / 9;
//      double g = x * x * x * x + y * y * y * y - x * x * y * y;
//      double g = (x * x * x * x * (13 * x * x - 15 * y * y) + y * y * y * y * (13 * y * y - 15 * x * x)) / 360;
//      double t = std::sqrt(std::max(std::fabs(x), std::fabs(y)));
      double g = r * r / 9;
      double t = std::min(x * x * x * x + y * y * y * y, 1.0);
      double val = t * f + (1 - t) * g;
#elif 1
      double val = r * r * r / 9;
      if (std::fabs(x) < 1 && std::fabs(y) < 1)
        val *= 1 - std::exp(-8 * (x * x * x * x + y * y * y * y - x * x * y * y));
//        val -= std::exp(-16 * (x * x * x * x + y * y * y * y)) / 64;
#else
      double val = r;
#endif
      set(u, n, i, j, val);
    }
  }

  return u;
}

// initialize from file
double*
init(size_t m, int g, const char* path)
{
  const size_t n = g + m + g; // number of grid cells including ghost layers

  FILE* file = fopen(path, "rb");
  if (!file) {
    fprintf(stderr, "cannot open file\n");
    return 0;
  }

  double* u = new double[n * n];
  if (fread(u, sizeof(*u), n * n, file) != n * n) {
    fprintf(stderr, "cannot read file\n");
    delete[] u;
    u = 0;
  }

  fclose(file);

  return u;
}

int main(int argc, char* argv[])
{
  const size_t m = SIZE; // number of grid cells per dimension
  const int g = GHOST; // number of ghost layers
  const size_t n = g + m + g; // number of grid cells including ghost layers
  const double h = 2. / m; // grid spacing for domain [-1, 1]^3
  size_t iters = ITERS; // number of iterations
  size_t diter = 64; // interval between frames

  // initialize u
  double* u = (argc == 2 ? init(m, g, argv[1]) : init(m, g));
  if (!u)
    return EXIT_FAILURE;

  // compress and decompress u
  double time = now();
  store(u, n);
  fprintf(stderr, "time = %.3f ms\n", 1000 * (now() - time));

#if JACOBI
  // Jacobi: keep intermediate state, v
  double* v = new double[n * n];
  std::copy(u, u + n * n, v);
#else
  // Gauss-Seidel: update u immediately
  double* v = u;
#endif

  // output initial solution
  if (diter)
    write(u, m, g, 0);

  // solve using Jacobi or Gauss-Seidel iteration
  for (size_t iter = 0; iter < iters; iter++) {
    // advance solution one iteration
    double max = 0;
    for (size_t j = g; j < n - g; j++) {
      double y = h * (j - g - 0.5 * (m - 1));
      for (size_t i = g; i < n - g; i++) {
        double x = h * (i - g - 0.5 * (m - 1));
        double r = std::sqrt(x * x + y * y);
        // Poisson equation: div(grad(u)) = r
        // uxx ~ (u(i-1, j) - 2 u(i, j) + u(i+1, j)) / h^2
        // uyy ~ (u(i, j-1) - 2 u(i, j) + u(i, j+1)) / h^2
        // u(i, j) = (sx + sy - h^2 r) / 4
        double sx = get(u, n, i - 1, j) + get(u, n, i + 1, j);
        double sy = get(u, n, i, j - 1) + get(u, n, i, j + 1);
#if 1
        // Poisson equation: div(grad(u)) = r
        double val = (sx + sy - h * h * r) / 4;
#else
        // Laplace equation: div(grad(u)) = 0
        double val = (sx + sy) / 4;
#endif
        // compute maximum deviation from previous solution
        max = std::max(max, std::fabs(get(u, n, i, j) - val));
        set(v, n, i, j, val);
      }
    }
    fprintf(stderr, "%zu %g\n", iters - iter, max);

    // compress v and decompress to u
    store(u, n, v);

    // output solution
    if (diter && !((iter + 1) % diter))
      write(u, m, g, (iter + 1) / diter);
  }

  fprintf(stderr, "total time = %.3f m\n", (now() - time) / 60);

  // output checkpoint
  fwrite(u, sizeof(*u), n * n, stdout);

  // clean up
  delete[] u;
  if (v != u)
    delete[] v;

  return 0;
}
