#include <iostream>
#if defined(__x86_64__) || defined(_M_X64)
#include <x86intrin.h>
#endif

#if defined(__x86_64__) || defined(_M_X64)
uint64_t inline rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#elif defined(__GNUC__) && defined(__riscv)  && !defined(HAVE_TICK_COUNTER)
static __inline__ uint64_t rdtsc(void)
{
        uint64_t cycles;
        __asm__ volatile("rdcycle %0" : "=r"(cycles));
        return cycles;
}
#endif

int main(int argc, char* argv[])
{
    uint64_t c1 = rdtsc();
    printf("start cycle is %llu\n", c1);
    
    for (int i=0; i<10; i++)
        printf("%d ", i*5);
    printf("\n");

    uint64_t c2 = rdtsc();
    printf("end cycle is %llu\n", c2);
    printf("cycle difference is %llu\n", c2-c1);

}