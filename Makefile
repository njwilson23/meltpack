TARGETS = cross.o eval.o fft2d.o fitreg.o gcorr.o gnorm.o kvert.o sums.o esterr.o

all: $(TARGETS)

%.o: src/%.f
	gfortran -c -fPIC $^

clean:
	find . -name "*.o" -print -delete

.PHONY: all clean
