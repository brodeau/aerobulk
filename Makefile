# AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com), S. Falahat (sd.falahat@gmail.com)
# https://github.com/brodeau/aerobulk/

include make.macro

All: lib/libaerobulk.a bin/test_aerobulk.x bin/example_call_aerobulk.x bin/test_phymbl.x bin/cx_vs_wind_test.x bin/test_ice.x bin/test_aerobulk_ice.x 

sk: bin/test_aerobulk_buoy_series_skin.x

CPP: lib/libaerobulk_cxx.a bin/example_call_aerobulk_cxx.x


# bin/test_coef_n10.x
# bin/test_coef_no98.x

#-L$(DIR_FORT_LIB) $(LNK_FORT_LIB)

LIB = -L./lib -laerobulk 

LIB_SRC = src/mod_const.f90 \
	  src/mod_phymbl.f90 \
          src/mod_cs_coare3p0.f90 \
          src/mod_skin_coare.f90 \
          src/mod_skin_ecmwf.f90 \
	  src/mod_blk_coare3p0.f90 \
	  src/mod_blk_coare3p6.f90 \
          src/mod_blk_ncar.f90 \
	  src/mod_blk_ecmwf.f90 \
	  src/mod_blk_ice_an05.f90 \
	  src/mod_blk_ice_lu12.f90 \
	  src/mod_blk_ice_lu15.f90 \
	  src/mod_blk_neutral_10m.f90 \
          src/mod_aerobulk_compute.f90 \
          src/mod_aerobulk.f90

LIB_OBJ = $(LIB_SRC:.f90=.o)

#LIB_OBO = $(LIB_SRC:.f90=.o)
#LIB_OBJ = $(patsubst src%,obj%,$(LIB_OBO))
#          src/mod_wl_coare3p6.f90 \
#          src/mod_wl_ecmwf.f90 \


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

bin/test_coef_n10.x: src/test_coef_n10.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_coef_n10.f90 -o bin/test_coef_n10.x $(LIB)

bin/test_coef_no98.x: src/test_coef_no98.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_coef_no98.f90 -o bin/test_coef_no98.x $(LIB)

bin/test_phymbl.x: src/test_phymbl.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_phymbl.f90 -o bin/test_phymbl.x $(LIB)

bin/cx_vs_wind_test.x: src/cx_vs_wind_test.f90 lib/libaerobulk.a
	@mkdir -p bin dat
	$(FC) $(FF) src/cx_vs_wind_test.f90 -o bin/cx_vs_wind_test.x $(LIB)

bin/test_ice.x: src/test_ice.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_ice.f90 -o bin/test_ice.x $(LIB)

bin/test_aerobulk_ice.x: src/test_aerobulk_ice.f90 lib/libaerobulk.a
	@mkdir -p bin
	$(FC) $(FF) src/test_aerobulk_ice.f90 -o bin/test_aerobulk_ice.x $(LIB)





bin/example_call_aerobulk_cxx.x: src/example_call_aerobulk.cpp lib/libaerobulk.a lib/libaerobulk_cxx.a
	@mkdir -p bin dat
	$(CXX) $(CXXFLAGS) src/example_call_aerobulk.cpp -o bin/example_call_aerobulk_cxx.x $(LIB_CXX) $(LIB)

bin/test_aerobulk_buoy_series_skin.x: src/test_aerobulk_buoy_series_skin.f90 lib/libaerobulk.a mod/io_ezcdf.mod
	@mkdir -p bin
	$(FC) $(FF) src/io_ezcdf.o src/test_aerobulk_buoy_series_skin.f90 -o bin/test_aerobulk_buoy_series_skin.x $(LIB) -L$(NETCDF_DIR)/lib $(L_NCDF)


mod/io_ezcdf.mod: src/io_ezcdf.f90
	@mkdir -p mod
	$(FC) $(FF) -I$(NETCDF_DIR)/include -c src/io_ezcdf.f90 -o src/io_ezcdf.o



.f90.o: $(LIB_SRC) $(LIB_SRC_CXX)
	@mkdir -p mod
	$(FC) -c $(FF) $< -o $*.o

.cpp.o: $(LIB_SRC_CXX)
	@mkdir -p mod
	$(CXX) -c $(CXXFLAGS) $< -o $*.o

clean:
	rm -rf mod bin lib src/*.o *~ \#* dat *.svg *.png *.eps *.gp *.out *.nc

#/*.x  *.log *~ out  mod/* lib/* *.nc tmp.* \#* *.info  config.dat *.tmp


