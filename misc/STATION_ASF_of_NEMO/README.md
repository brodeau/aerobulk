
## WARNING: TOTALLY-ALPHA-STUFF / DOCUMENT IN THE PROCESS OF BEING WRITEN!

# *Station Air-Sea Fluxes* demonstration case

## Objectives

```STATION_ASF``` is a demonstration case that mimics an in-situ station (buoy, platform) dedicated to the estimation of surface air-sea fluxes by means of the measurement of traditional meteorological surface parameters.

```STATION_ASF``` is based on the merging of the "single column" and the "standalone surface module" configurations of NEMO. In short, it coulb defined as "SAS meets C1D". As such, the spatial domain of ```STATION_ASF``` is punctual (1D, well actually 3 x 3 as in C1D).

```STATION_ASF``` is therefore a versatile tool, and extremely light in terms of computing requirements, to test the different bulk algorithms and cool-skin/warm-layer parameterization options included in NEMO.

As input ```STATION_ASF``` will require the traditional *bulk* sea surface parameters:

- sea surface temperature (SST) at $z_{SST}$ meters below the surface
- Surface current vector
- Sea surface salinity

as well as the usual surface atmospheric state:

- air temperature at $z_t$ meters above the surface
- air humidity  at $z_t$ meters above the surface (specific humidity or relative humidity or dew-point temperature)
- wind speed vector at $z_u$ meters above the surface
- Sea level atmospheric pressure (SLP)
- Downwelling solar radiation
- Downwelling longwave radiation



## Physical description

### Important namelist parameters speficic to STATION_ASF

* ```rn_dept1@namusr_def:``` depth (m) at which the prescribed SST is taken (i.e. depth of first T-point); important due to impact on warm-layer estimate, the deeper, the more pronounced!

* ```rn_lat1d,rn_lon1d@namc1d:``` fixed coordinates of the location of the station (buoy, platform, etc).

* ```namsbc_blk:``` to be filled carefully, just as for "C1D", the prescribed surface ATMOSPHERIC state (files) are time series of shape 3x3 in space

* ```namsbc_sas:``` to be filled carefully, just as for "C1D", the prescribed surface OCEAN state (files) are time series of shape 3x3 in space



## Input files to test STATION ASF

Three full years of processed hourly data from the PAPA station (buoy) can be downloaded here:
https://drive.google.com/file/d/1MxNvjhRHmMrL54y6RX7WIaM9-LGl--ZP/

These three files are everything you need to play with the set of namelists provided for this test-case.

- ```Station_PAPA_50N-145W_atm_hourly.nc```  → contains hourly surface atmospheric state
- ```Station_PAPA_50N-145W_precip_daily.nc``` → contains daily precipitation
- ```Station_PAPA_50N-145W_oce_hourly.nc``` → contains hourly sea surface state

For station PAPA (50.1 N, 144.9 W), air temperature and humidity are measured at 2.5 m, the wind speed at 4 m, and the SST at 1 m below the surface, hence the following namelist parameters are given:

- ```rn_dept1 =    1.  ``` (&namusr_def)
- ```rn_lat1d =  50.1 ``` (&namc1d)
- ```rn_lon1d = 215.1``` (&namc1d)
- ```rn_zqt   =   2.5``` (&namsbc_blk)
- ```rn_zu    =    4.``` (&namsbc_blk)



## Playing with STATION_ASF

First compile the test-case as follows (compile with xios-2.5 support → check your ARCH file):

```./makenemo -m <your_arch> -n STATION_ASF -j 4 -a STATION_ASF```

Then you can use the script ``launch_sasf.sh`` found in  ```EXPREF/``` to launch 3 simulations (one for each bulk parameterization available). You need to adapt the following variable to your environment in the script:

- ```NEMO_DIR``` : NEMO root directory where to fetch compiled STATION_ASF ```nemo.exe``` + setup (such as ```${NEMO_DIR}/tests/STATION_ASF```)

- ```WORK_DIR``` :  Directory where to run the simulation

- ```FORC_DIR```  Directory containing sea-surface + atmospheric forcings (get it there https://drive.google.com/file/d/1MxNvjhRHmMrL54y6RX7WIaM9-LGl--ZP/)

