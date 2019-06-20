ifeq ($(PREFIX),)
	PREFIX := build
endif


F90 = gfortran
#FFLAGS = -Wall -ffpe-trap=invalid,zero -O0 -fPIC -g -pg
FFLAGS = -Wall -ffpe-trap=invalid,zero -O3 -fPIC

all: $(addprefix build/, test_compton test_external_compton test_ec_blr test_ec_dust test_absorption liboptical_depth.a)
objects_compton = $(addprefix build/, const.o quadpack.o compton.o dilog.o)
objects_external = $(addprefix build/, gamma_avg.o zeroin.o external_compton.o)
objects_ssc = $(addprefix build/, const.o gamma_avg.o dilog.o ssc.o zeroin.o quadpack.o)

# Linking of the executables
build/liboptical_depth.a: $(addprefix build/, optical_depth.o const.o photoabsorption.o multidim_integrate.o quadpack.o optical_depth.mod const.mod photoabsorption.mod multidim_integrate.mod)
	ar cr $@ $^
build/test_compton: build/test_compton.o $(objects_compton) | build
	${F90} ${FFLAGS} -o $@ $^
build/test_external_compton: build/test_external_compton.o $(objects_compton)\
 $(objects_external) | build
	${F90} ${FFLAGS} -o $@ $^
build/test_ec_blr: build/test_ec_blr.o $(objects_compton) $(objects_external)\
 | build
	${F90} ${FFLAGS} -o $@ $^
build/test_ec_dust: build/test_ec_dust.o $(objects_compton) $(objects_external)\
 | build
	${F90} ${FFLAGS} -o $@ $^

# Add Module dependencies
build/optical_depth.o: $(addprefix build/, const.o photoabsorption.o)
build/test_compton.o: $(addprefix build/, compton.o)
build/test_external_compton.o: $(addprefix build/, external_compton.o)
build/test_ec_blr.o: $(addprefix build/, external_compton.o)
build/test_ec_dust.o: $(addprefix build/, external_compton.o)

build/compton.o: $(addprefix build/, const.o quadpack.o dilog.o)
build/gamma_avg.o: $(addprefix build/, zeroin.o compton.o)
build/external_compton.o: $(addprefix build/, compton.o const.o gamma_avg.o)
build/ssc.o: $(addprefix build/, const.o gamma_avg.o dilog.o)
build/photoabsorption.o: $(addprefix build/, const.o multidim_integrate.o quadpack.o)

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

# Stuff for testing
build/3C454_3: build/3C454_3.o $(objects_compton) $(objects_external)\
 $(objects_ssc) build/luminosity_distance.o build/photoabsorption.o build/multidim_integrate.o | build
	${F90} ${FFLAGS} -o $@ $^
build/4C71_07: build/4C71_07.o $(objects_compton) $(objects_external)\
 $(objects_ssc) build/luminosity_distance.o build/photoabsorption.o build/multidim_integrate.o | build
	${F90} ${FFLAGS} -o $@ $^
build/dermer: build/dermer.o $(objects_compton) $(objects_external)\
 $(objects_ssc) build/luminosity_distance.o build/photoabsorption.o build/multidim_integrate.o | build
	${F90} ${FFLAGS} -o $@ $^
build/dermer.o: $(addprefix build/, ssc.o external_compton.o photoabsorption.o  luminosity_distance.o)
build/W-Com.o: $(addprefix build/, ssc.o external_compton.o photoabsorption.o luminosity_distance.o)
build/4C71_07.o: $(addprefix build/,ssc.o external_compton.o photoabsorption.o  luminosity_distance.o)

build/test_absorption.o: $(addprefix build/, photoabsorption.o const.o)
build/test_absorption: $(addprefix build/, test_absorption.o photoabsorption.o multidim_integrate.o const.o quadpack.o)
	${F90} ${FFLAGS} -o $@ $^
build/3C84.o: $(addprefix build/, photoabsorption.o const.o)
build/3C84: $(addprefix build/, 3C84.o photoabsorption.o multidim_integrate.o const.o quadpack.o)
	${F90} ${FFLAGS} -o $@ $^
