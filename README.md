

![Aerobulk Logo](https://brodeau.github.io/images/projects/aerobulk_logo_s.svg)

<!-- Online documentation: https://brodeau.github.io/aerobulk/ -->

**AeroBulk** is a FORTRAN90-based library and suite of tools that feature the most accurate parameterizations available to estimate turbulent air-sea fluxes by means of the traditional **aerodynamic bulk formulae**.

These turbulent fluxes, namely, wind stress, evaporation (latent heat flux) and
sensible heat flux, are estimated using the sea surface temperature (bulk or
skin), and the near-surface atmospheric surface state: wind speed, air
temperature and humidity.

![Bulk Formula](https://brodeau.github.io/images/projects/bulk.svg)


The following figure provides a schematic view on the way turbulent fluxes are computed in AeroBulk:

![Aerobulk Approach](https://brodeau.github.io/images/projects/fig_bulk_model_ftp.svg)

&nbsp;

In AeroBulk, 5 *state of the art* bulk algorithms are available to compute the drag, sensible heat and moisture transfer coefficients (namely C<sub>D</sub>, C<sub>H</sub> and C<sub>E</sub>) used in the bulk formula:

*   COARE v3.0 ([Fairall *et al.*, 2003](http://dx.doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2))
*   COARE v3.6 (Fairall *et al.*, 2018 + [Edson *et al.*, 2013](http://dx.doi.org/10.1175/jpo-d-12-0173.1))
*   ECMWF ([IFS (Cy40) documentation](https://software.ecmwf.int/wiki/display/IFS/CY40R1+Official+IFS+Documentation))
*   ANDREAS ([Andreas *et al.*, 2015](https://dx.doi.org/10.1002/qj.2424))
*   NCAR (Large & Yeager 2004, [2009](http://dx.doi.org/10.1007/s00382-008-0441-3))

In the COARE and ECMWF algorithms, a cool-skin/warm layer scheme is included and can be activated if the input sea-surface temperature is the bulk SST (usually measured a few tenths of meters below the surface). Activation of these cool-skin/warm layer schemes requires the surface downwelling shortwave and longwave radiative flux components to be provided. The NCAR algorithm is supposed to be used with the bulk SST and does not feature a cool-skin/warm layer scheme.

Beside bulk algorithms, AeroBulk also provides a variety of functions to accurately estimate relevant atmospheric variable such as density of air, different expressions of the humidity of air, viscosity of air, specific humidity at saturation, Obukhov length, wind gustiness, etc...

The focus in AeroBulk is readability, efficiency, and portability towards modern ocean & atmosphere GCMs (Fortran 90, set of modules and a library).

&nbsp;

![Aerobulk Approach](https://brodeau.github.io/images/projects/Comparaison_Psi.svg)
*Example of a set of figures generated with one of AeroBulk diagnostic scripts: comparison of the stability correction profiles Psi(zeta) as used in 4 different bulk algorithms.*

&nbsp;


# **> Giving AeroBulk a first try in interactive "toy mode"**

Check that ```bin/aerobulk_toy.x``` has been compiled and execute it:

     ./bin/aerobulk_toy.x

You will be guided and interactively prompted for different sea-surface and ABL related parameters, such as SST, air temperature and humidity, wind speed, etc.  Then ```aerobulk_toy.x``` will compute all the turbulent air-sea fluxes with all the algorithms available (as well as third-party diagnostics of the ABL), and will print it in the form of a summary table.

 List of command line options for ```aerobulk_toy.x```:

    -p   => Ask for sea-level pressure, otherwise assume 1010 hPa
    
    -r   => Ask for relative humidity rather than specific humidity
    
    -S   => Use the Cool Skin Warm Layer parameterization to compute
            and use the skin temperature instead of the bulk SST
            only in COARE and ECMWF
    
    -N   => Force neutral stability in surface atmospheric layer
            -> will not ask for air temp. at zt, instead will only
            ask for relative humidity at zt and will force the air
            temperature at zt that yields a neutral-stability!
    
    -h   => Show this message


Example of an output obtained with the following setup:

 * _z<sub>u_ = 10 m
 * _z<sub>t_ = 2 m
 * _SST_ = 22&deg;C
 * _t<sub>zt_ = 20&deg;C
 * _q<sub>zt_ = 12 g/kg
 * _U<sub>zu_ = 5 m/s

<!-- $$  SST = 22^{\circ}C \\ \theta_z = 20^{\circ}C $$ -->

     ==============================================================================================
        Algorithm:      coare3p0  |  coare3p6  |    ncar    |    ecmwf   | andreas
     ==============================================================================================
           C_D     =    1.1954       1.0775       1.2038       1.2862       1.0167      [10^-3]
           C_E     =    1.3345       1.3729       1.3618       1.3143       1.1565      [10^-3]
           C_H     =    1.3345       1.3729       1.2776       1.2635       1.1103      [10^-3]
    
           z_0     =   4.40936E-05  2.19285E-05  4.49880E-05  6.98835E-05  1.56119E-05  [m]
           u*      =   0.17578      0.16672      0.17348      0.18192       0.1594      [m/s]
           L       =   -20.383      -16.919      -20.494      -24.029      -18.558      [m]
           Ri_bulk =  -3.78706E-02 -3.79537E-02 -3.90686E-02 -3.79799E-02 -3.87826E-02  [-]
    
                      *** Neutral stability: **
           UN10    =    5.4192       5.4311       5.3396       5.3992       5.3289      [m/s]
           C_D_N   =    1.0521      0.94234       1.0555       1.1353       0.8950      [10^-3]
           C_E_N   =    1.1077       1.1119       1.1241       1.1064       0.9600      [10^-3]
           C_H_N   =    1.1077       1.1119       1.0624       1.0680       0.9260      [10^-3]
    
     Equ. Charn p. =   1.10033E-02  4.23674E-03  1.15475E-02  1.80029E-02  2.02298E-03
    
      Wind stress  =    36.198       32.596       35.849       38.860       30.275      [mN/m^2]
      Evaporation  =    3.1059       3.1925       3.1215       3.0538       2.6270      [mm/day]
         QL        =   -88.031      -90.486      -88.473      -86.555      -74.459      [W/m^2]
         QH        =   -17.833      -18.330      -16.730      -16.805      -14.441      [W/m^2]


&nbsp;

# **> Computing turbulent fluxes with AeroBulk**

AeroBulk can also directly compute the 3 turbulent fluxes with the routine ```aerobulk_model()``` of module ```mod_aerobulk``` (mod\_aerobulk.f90):

       PROGRAM TEST_FLUX
           USE mod_aerobulk
           ...
           CALL AEROBULK_MODEL( calgo, zt, zu, sst, t_zt, q_zt, U_zu, V_zu, SLP, &
           &                    Qe, Qh, Tau_x, Tau_y                             &
           &                   [, Niter=N, rad_sw=Rsw, rad_lw=Rlw] )
           ...
       END PROGRAM TEST_FLUX

INPUT ARGUMENTS:

*   ```calgo``` : (String) algorithm to use (coare3p0/coare3p6/ncar/ecmwf)
*   ```zt``` : (Sc,real) height for temperature and spec. hum. of air [m]
*   ```zu``` : (Sc,real) height for wind speed (generally 10m) [m]
*   ```sst``` : (2D,real) SST [K]
*   ```t_zt``` : (2D,real) potential air temperature at zt [K]
*   ```q_zt``` : (2D,real) specific humidity of air at zt [kg/kg]
*   ```U_zu``` : (2D,real) zonal scalar wind speed at 10m [m/s]
*   ```V_zu``` : (2D,real) meridional scalar wind speed at 10m [m/s]
*   ```SLP``` : (2D,real) sea-level pressure [Pa]

[ OPTIONAL INPUT ARGUMENT: ]

*   ```Niter``` : (Sc,int) number of iterations (default is 4)
*   ```rad_sw``` : (2D,real) downw. shortwave rad. at surface (>0) [W/m^2]
*   ```rad_lw``` : (2D,real)downw. longwave rad. at surface (>0) [W/m^2]

(The presence of ```rad_sw``` and ```rad_sw``` triggers the use of the Cool-Skin Warm-Layer parameterization with COARE and ECMWF algorithms)

OUTPUT ARGUMENTS:

*   ```Qe``` : (2D,real) latent heat flux [W/m^2]
*   ```Qh``` : (2D,real) sensible heat flux [W/m^2]
*   ```Tau_x``` : (2D,real) zonal wind stress [N/m^2]
*   ```Tau_y``` : (2D,real) meridional wind stress [N/m^2]

Example of a call, using COARE 3.0 algorithm with cool-skin warm-layer parameterization and 10 iterations:

           CALL AEROBULK_MODEL( 'coare_3p0', 2., 10., sst, t_zt, q_zt, U_zu, V_zu, SLP, &
           &                    Qe, Qh, Tau_x, Tau_y,                                   &
           &                    Niter=10, rad_sw=Rsw, rad_lw=Rlw )

&nbsp;

# **> Computing transfer coefficients with AeroBulk**

In AeroBulk, 3 different routines are available to compute the bulk transfer (_a.k.a_ exchange) coefficients C<sub>D</sub>, C<sub>H</sub> and C<sub>E</sub>. Beside computing the transfer coefficients, these routines adjust air temperature and humidity from height _z<sub>t</sub>_ to the reference height (wind) _z<sub>u</sub>_. They also return the bulk wind speed, which is the scalar wind speed at height _z<sub>u</sub>_ with the potential inclusion of a gustiness contribution (in calm and unstable conditions).


**> TURB_COARE_3p0, transfer coefficients with COARE 3.0 (replace "3p0" by "3p6" to use COARE 3.6)**

Use ```turb_coare_3p0()``` of module ```mod_blk_coare_3p0``` (mod\_blk\_coare_3p0.f90).
Example of a call:

              PROGRAM TEST_COEFF
                  USE mod_const
                  USE mod_blk_coare_3p0
                  ...
                  jpi = Ni ! x-shape of the 2D domain
                  jpj = Nj ! y-shape of the 2D domain
                  ...
                  CALL TURB_COARE_3P0( zt, zu, T_s, t_zt, q_s, q_zt, U_zu,  &
                  &                Cd, Ch, Ce, t_zu, q_zu, U_blk            &
                  &                [ , Qsw=Rsw, rad_lw=Rlw, slp=P ]         &
                  &                [ , xz0=z0, xu_star=u_s, xL=L ] )
                  ...
              END PROGRAM TEST_COEFF

**INPUT ARGUMENTS:**

*   ```zt``` : (Sc,real) height for air temperature and humidity [m]
*   ```zu``` : (Sc,real) height for wind speed (generally 10m) [m]
*   ```t_zt``` : (2D,real) potential air temperature at zt [K]
*   ```q_zt``` : (2D,real) air spec. humidity of at zt [kg/kg]
*   ```U_zu``` : (2D,real) scalar wind speed at zu [m/s]

**INPUT and OUTPUT ARGUMENTS:**

*   ```T_s``` : (2D,real) surface temperature [K]
    *   input: bulk SST
    *   output: skin temperature or SST (unchanged)
*   ```q_s``` : (2D,real) surface satur. spec. humidity at T_s [kg/kg]
    *   input: saturation at bulk SST (not needed if skin p. used)
    *   output: saturation at skin temp. or at SST (unchanged)

**[ OPTIONAL INPUT ARGUMENTS: ]**

*   ```Qsw``` : (2D,real) net shortw. rad. at surface (>0, after albedo!) [W/m^2]
*   ```rad_lw``` : (2D,real) downw. longw. rad. at surface (>0) [W/m^2]
*   ```slp``` : (2D,real) sea-level pressure [Pa]

(The presence of these 3 optional input parameters triggers the use of the Cool-Skin Warm-Layer parameterization)

**OUTPUT ARGUMENTS:**

*   ```Cd``` : (2D,real) drag coefficient
*   ```Ch``` : (2D,real) sensible heat transfer coefficient
*   ```Ce``` : (2D,real) moisture transfer (evaporation) coefficient
*   ```t_zu``` : (2D,real) air pot. temperature adjusted at zu [K]
*   ```q_zu``` : (2D,real) air spec. humidity adjusted at zu [kg/kg]
*   ```Ublk``` : (2D,real) bulk wind speed at 10m [m/s]

**[ OPTIONAL OUTPUT ARGUMENTS: ]**

*   ```z0``` : (2D,real) roughness length of the sea surface [m]
*   ```u_s``` : (2D,real) friction velocity [m/s]
*   ```L``` : (2D,real) Obukhov length [m]

**> Some Examples**

Using COARE 3.0 without the cool-skin warm-layer parameterization, with air temperature and humidity provided at 2m and wind at 10m:

              PROGRAM TEST_COEFF
                  USE mod_const
                  USE mod_blk_coare_3p0
                  ...
                  jpi = Ni ! x-shape of the 2D domain
                  jpj = Nj ! y-shape of the 2D domain
                  ...
                  CALL TURB_COARE_3P0( '3.0', 2., 10., Ts, t2, qs, q2, U10, &
                  &                Cd, Ch, Ce, t10, q10, U_blk )
                  ...
              END PROGRAM TEST_COEFF

In this case, Ts and qs, the surface temperature and saturation specific humidity won't be modified. The relevant value of qs must be provided as input.

Now the same but using the cool-skin warm-layer parameterization:

              PROGRAM TEST_COEFF
                  USE mod_const
                  USE mod_blk_coare_3p0
                  ...
                  jpi = Ni ! x-shape of the 2D domain
                  jpj = Nj ! y-shape of the 2D domain
                  ...
                  CALL TURB_COARE_3P0( '3.0', 2., 10., Ts, t2, qs, q2, U10, &
                  &                Cd, Ch, Ce, t10, q10, U_blk,         &
                  &                Qsw=Rsw, rad_lw=Rlw, slp=MSL )
                  ...
              END PROGRAM TEST_COEFF

Here, Ts is the bulk SST as input and will become the skin temperature as output! qs is irrelevant as input and is the saturation specific humidity at temperature Ts as output!


&nbsp;

# **> Computing atmospheric state variables with AeroBulk**

A selection of useful functions to estimate some atmospheric state variables in the marine boundary layer are available in the module ```mod_phymbl``` (mod\_phymbl.f90).

Example for computing SSQ of Eq.(1) out of the SST and the SLP:

              PROGRAM TEST_PHYMBL
                  USE mod_const
                  USE mod_phymbl
                  ...
    
                  SSQ(:,:) = 0.98*q_sat(SST, SLP)
                  ...
              END PROGRAM TEST_PHYMBL




&nbsp;

**> Acknowledging AeroBulk**

_To acknowledge/reference AeroBulk in your scientific work, please cite:_
Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically significant effects of some approximations in the bulk parameterizations of turbulent air-sea fluxes. _J. Phys. Oceanogr._, **47 (1)**, 5–28, 10.1175/JPO-D-16-0169.1. [ [DOI](http://dx.doi.org/10.1175/JPO-D-16-0169.1) ]
