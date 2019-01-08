ifeq ($(PREFIX),)
	PREFIX := build
endif


F90 = gfortran
FFLAGS = -Wall -ffpe-trap=invalid,zero -O3 -g -fPIC

all: $(addprefix build/, test_compton test_external_compton)
objects_compton = $(addprefix build/, const.o quadpack.o compton.o dilog.o)
objects_external = $(addprefix build/, gamma_avg.o zeroin.o external_compton.o)

# Linking of the executables
build/test_compton: build/test_compton.o $(objects_compton) | build
	${F90} ${FFLAGS} -o $@ $^
build/test_external_compton: build/test_external_compton.o $(objects_compton)\
 $(objects_external) | build
	${F90} ${FFLAGS} -o $@ $^

# Add Module dependencies
build/test_compton.o: $(addprefix build/, compton.o)
build/test_external_compton.o: $(addprefix build/, external_compton.o)

build/compton.o: $(addprefix build/, const.o quadpack.o dilog.o)
build/gamma_avg.o: $(addprefix build/, zeroin.o compton.o)
build/external_compton.o: $(addprefix build/, compton.o const.o gamma_avg.o)

build/%.o build/%.mod: %.f | build
	${F90} ${FFLAGS} -J build -o build/$*.o -c $<

# pattern rule for all fortran source files
build/%.o build/%.mod: %.f90 | build
	${F90} ${FFLAGS} -J build -o build/$*.o -c $<

build:
	mkdir -p $@

clean:
	rm -rf build

.PHONY: all clean install
