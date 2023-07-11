# Instructions for running diffusion w compression software cache
## Compiling
```
// From the repo root compile zfp library
make
// Go to example folder to compilethe diffusion benchmark
cd examples/ && make ../bin/diffusion
```

## Run command
```
./bin/diffusion -n 128 128 -t 10 -r 64 >> 128_128_10_64.out
./bin/diffusion -n 128 128 -t 10 -r 16 >> 128_128_10_16.out
```