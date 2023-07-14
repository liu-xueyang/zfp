#DEFS = -DREAL=__fp16
#DEFS = -DREAL=float
#DEFS = -DREAL=double
#DEFS = "-DREAL=NumRep<unsigned short, ExpGolombRice<2>, MapLinear>"
#DEFS = "-DREAL=NumRep<unsigned int, ExpGolombRice<2>, MapLinear>"

#ZFP_RATE = 5.5
#ZFP_RATE = 6
#ZFP_RATE = 8
#ZFP_RATE = 10
#ZFP_RATE = 12
ZFP_RATE = 16
DEFS = -DZFP_RATE=$(ZFP_RATE)

#DEFS = -DBLAZ

# output gradient norm rather than Laplacian
#DEFS += -DGRADIENT

# output function rather than Laplacian
#DEFS += -DFUNCTION

DEFS += -DSIZE=1024
DEFS += -DITERS=64
#DEFS += -DSIZE=4096
#DEFS += -DITERS=0

all: poisson

poisson: poisson.cpp
	$(CROSS_COMPILE)g++ -O3 $< -o $@ $(DEFS) -I../include -L../lib -lzfp

exact: exact.cpp
	$(CROSS_COMPILE)g++ -O3 $< -o $@

clean:
	rm -f poisson exact
