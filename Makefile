TARGETS = src/cross.o src/eval.o src/fft2d.o src/fitreg.o src/gcorr.o src/gnorm.o src/kvert.o src/sums.o

all: $(TARGETS)

%.o: %.f
	gfortran -c -fPIC $^

clean:
	find src/ -name "*.o" -print -delete

.PHONY: all clean
