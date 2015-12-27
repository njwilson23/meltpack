TARGETS = cross.o eval.o fft2d.o fitreg.o gcorr.o gnorm.o kvert.o sums.o

all: $(TARGETS)

%.o: %.f
	gfortran -c $^
	cp $@ MELTPACK/$@

clean:
	find . -name "*.o" -print -delete

.PHONY: all clean
