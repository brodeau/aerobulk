# AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com), S. Falahat (sd.falahat@gmail.com)
# https://sourceforge.net/p/aerobulk

include make.macro

All: lib/libaerobulk.a bin/test_aerobulk.x bin/test_skin_corr.x bin/test_coef_skin.x \
	bin/example_call_aerobulk.x bin/test_thermo.x bin/cx_vs_wind_test.x lib/libaerobulk_cxx.a bin/example_call_aerobulk_cxx.x

# bin/test_coef_n10.x
# bin/test_coef_no98.x

LIB = -L./lib -laerobulk

LIB_SRC = src/mod_const.f90 \
	  src/mod_thermo.f90 \
	  src/mod_blk_neutral_10m.f90 \
	  src/mod_blk_coare.f90 \
          src/mod_blk_ncar.f90 \
	  src/mod_blk_ecmwf.f90 \
          src/mod_aerobulk_compute.f90 \
          src/mod_aerobulk.f90

LIB_OBJ = $(LIB_SRC:.f90=.o)

LIB_CXX = -L./lib -laerobulk_cxx

LIB_SRC_CXX = src/aerobulk.cpp \
		  src/aerobulk_cxx.f90

LIB_OBJ_CXX = src/aerobulk.o \
		  src/aerobulk_cxx.o


CXXFLAGS += -I./include

.SUFFIXES: 
.SUFFIXES: .f90 .o .cpp

lib/libaerobulk.a: $(LIB_OBJ)
	@echo ""
	@mkdir -p lib
	ar -rv lib/libaerobulk.a  $(LIB_OBJ)
	ranlib lib/libaerobulk.a
	@echo ""

lib/libaerobulk_cxx.a: $(LIB_OBJ) $(LIB_OBJ_CXX)
	@echo ""
	@mkdir -p lib
	ar -rv lib/libaerobulk_cxx.a $(LIB_OBJ_CXX)
	ranlib lib/libaerobulk_cxx.a
	@echo ""

bin/test_aerobulk.x: src/test_aerobulk.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_aerobulk.f90 -o bin/test_aerobulk.x $(LIB)

bin/example_call_aerobulk.x: src/example_call_aerobulk.f90 lib/libaerobulk.a #
	@mkdir -p bin
	$(FC) $(FF) src/example_call_aerobulk.f90 -o bin/example_call_aerobulk.x $(LIB)

bin/test_skin_corr.x: src/test_skin_corr.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_skin_corr.f90 -o bin/test_skin_corr.x $(LIB)

bin/test_coef_n10.x: src/test_coef_n10.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_coef_n10.f90 -o bin/test_coef_n10.x $(LIB)

bin/test_coef_skin.x: src/test_coef_skin.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_coef_skin.f90 -o bin/test_coef_skin.x $(LIB)

bin/test_coef_no98.x: src/test_coef_no98.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_coef_no98.f90 -o bin/test_coef_no98.x $(LIB)

bin/test_thermo.x: src/test_thermo.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_thermo.f90 -o bin/test_thermo.x $(LIB)

bin/cx_vs_wind_test.x: src/cx_vs_wind_test.f90 lib/libaerobulk.a
	@mkdir -p bin dat
	$(FC) $(FF) src/cx_vs_wind_test.f90 -o bin/cx_vs_wind_test.x $(LIB)

bin/example_call_aerobulk_cxx.x: src/example_call_aerobulk.cpp lib/libaerobulk.a lib/libaerobulk_cxx.a
	@mkdir -p bin dat
	$(CPP) $(CXXFLAGS) src/example_call_aerobulk.cpp -o bin/example_call_aerobulk_cxx.x $(LIB_CXX) $(LIB)

.f90.o: $(LIB_SRC) $(LIB_SRC_CXX)
	@mkdir -p mod
	$(FC) -c $(FF) $< -o $*.o

.cpp.o: $(LIB_SRC_CXX)
	@mkdir -p mod
	$(CPP) -c $(CXXFLAGS) $< -o $*.o

clean:
	rm -rf mod bin lib src/*.o *~ \#* dat *.svg *.png

#/*.x  *.log *~ out  mod/* lib/* *.nc tmp.* \#* *.info  config.dat *.tmp


