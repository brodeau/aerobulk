<center><embed type="image/png" src="images/aerobulk_logo.png" style="width:500px"></center>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 20px;width: 960px;" align="justify">**AeroBulk** is a package/library that gathers state-of-the-art aerodynamic bulk formulae algorithms used to estimate turbulent air-sea fluxes in an efficient and unified way. These turbulent fluxes are wind stress, evaporation (latent heat flux) and sensible heat flux, they are needed as part of the surface boundary conditions of OGCMs, AGCMs and in the coupling interface of Earth Systems.  

AeroBulk relies on bulk formulae to compute turbulent air-sea fluxes from the sea surface temperature, wind speed, and air temperature and specific humidity. In AeroBulk, 4 state-of-the-art algorithms are available to compute the drag, sensible heat and moisture transfer coefficients (C<sub>D</sub>, C<sub>H</sub> and C<sub>E</sub>) used in the bulk formulaes:

*   COARE v3.0 ([Fairall _et al._ 2003](http://dx.doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2))
*   COARE v3.5 ([Edson _et al._ 2013](http://dx.doi.org/10.1175/jpo-d-12-0173.1)<a>)</a>
*   <a>ECMWF (</a>[IFS (Cy40) documentation](https://software.ecmwf.int/wiki/display/IFS/CY40R1+Official+IFS+Documentation))
*   NCAR (Large & Yeager 2004, [2009](http://dx.doi.org/10.1007/s00382-008-0441-3))

In the COARE and ECMWF algorithms, a cool-skin/warm layer scheme is included and can be activated if the input sea-surface temperature is the bulk SST (measured a few tenths of meters below the surface). Activation of these cool-skin/warm layer schemes requires the surface downwelling shortwave and longwave radiative flux components to be provided. The NCAR algorithm is to be used only with the bulk SST.  

Beside bulk algorithms AeroBulk also provides a variety of functions to accurately estimate relevant atmospheric variable such as density of air, different expressions of the humidity of air, viscosity of air, specific humidity at saturation, Monin-Obukhov length, wind gustiness, etc...  

The focus in AeroBulk is readability, efficiency and portability towards either modern GCMs (Fortran 90, set of modules and a library).</div>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 24px;width: 960px;" align="left">**> Obtaining AeroBulk**</div>

<div style="margin: 10px 40px 0px 60px; font-family: Courier New; font-size: 20px;width: 960px;" align="justify">The AeroBulk project is hosted by SourceForge.net and the source code is available here: [http://sourceforge.net/projects/aerobulk/](http://sourceforge.net/projects/aerobulk/)  

Alternatively you can directly checkout the code via subversion:

      svn co svn://svn.code.sf.net/p/aerobulk/code/trunk aerobulk

</div>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 24px;width: 960px;" align="left">**> Computing transfer coefficients with AeroBulk**</div>

<div style="margin: 10px 40px 0px 60px; font-family: Courier New; font-size: 20px;width: 960px;" align="justify">In AeroBulk, 3 different routines are available to compute the bulk transfer (_a.k.a_ exchange) coefficients C<sub>D</sub>, C<sub>H</sub> and C<sub>E</sub>. Beside computing the transfer coefficients, these routines adjust air temperature and humidity from height _z<sub>t</sub>_ to the reference height (wind) _z<sub>u</sub>_. They also return the bulk wind speed, which is the scalar wind speed at height _z<sub>u</sub>_ with the potential inclusion of a gustiness contribution (in calm and unstable conditions).

*   _[turb_coare()](src/turb_coare.html)_ of module **mod_blk_coare** (mod_blk_coare.f90)
*   _[turb_ecmwf()](src/turb_ecmwf.html)_ of module **mod_blk_ecmwf** (mod_blk_ecmwf.f90)
*   _[turb_ncar()](src/turb_ncar.html)_ of module **mod_blk_ncar** (mod_blk_ncar.f90)

</div>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 24px;width: 960px;" align="left">**> Computing turbulent fluxes with AeroBulk**</div>

<div style="margin: 10px 40px 0px 40px; font-family: Courier; font-size: 20px;width: 960px;" align="justify">AeroBulk can also directly compute the 3 turbulent fluxes with the routine _aerobulk_model()_ of module **mod_aerobulk** (mod_aerobulk.f90):

       PROGRAM TEST_FLUX
           USE mod_aerobulk
           ...
           CALL AEROBULK_MODEL( calgo, zt, zu, sst, t_zt, q_zt, U_zu, V_zu, SLP, &
           &                    Qe, Qh, Tau_x, Tau_y                             &
           &                   [, Niter=N, rad_sw=Rsw, rad_lw=Rlw] )
           ...
       END PROGRAM TEST_FLUX

INPUT ARGUMENTS:

*   calgo: (String) algorithm to use (coare/coare35/ncar/ecmwf)
*   zt : (Sc,real) height for temperature and spec. hum. of air [m]
*   zu : (Sc,real) height for wind speed (generally 10m) [m]
*   sst : (2D,real) SST [K]
*   t_zt : (2D,real) potential air temperature at zt [K]
*   q_zt : (2D,real) specific humidity of air at zt [kg/kg]
*   U_zu : (2D,real) zonal scalar wind speed at 10m [m/s]
*   V_zu : (2D,real) meridional scalar wind speed at 10m [m/s]
*   SLP : (2D,real) sea-level pressure [Pa]

[ OPTIONAL INPUT ARGUMENT: ]

*   Niter: (Sc,int) number of iterations (default is 4)
*   rad_sw: (2D,real) downw. shortwave rad. at surface (>0) [W/m^2]
*   rad_lw: (2D,real)downw. longwave rad. at surface (>0) [W/m^2]

<div style="margin: -18px 40px 0px 20px; font-family: Courier; font-size: 18x;width: 860px;" align="left">(The presence of rad_sw and rad_sw triggers the use of the Cool-Skin Warm-Layer parameterization with COARE* and ECMWF algorithms)</div>

OUTPUT ARGUMENTS:

*   Qe : (2D,real) latent heat flux [W/m^2]
*   Qh : (2D,real) sensible heat flux [W/m^2]
*   Tau_x : (2D,real) zonal wind stress [N/m^2]
*   Tau_y : (2D,real) meridional wind stress [N/m^2]

Example of a call, using COARE 3.0 algorithm with cool-skin warm-layer parameterization and 10 iterations:

           CALL AEROBULK_MODEL( 'coare', 2., 10., sst, t_zt, q_zt, U_zu, V_zu, SLP, &
           &                    Qe, Qh, Tau_x, Tau_y,                               &
           &                    Niter=10, rad_sw=Rsw, rad_lw=Rlw )

</div>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 24px;width: 960px;" align="left">**> Computing atmospheric state variables with AeroBulk**</div>

<div style="margin: 10px 40px 0px 60px; font-family: Courier New; font-size: 20px;width: 960px;" align="justify">A selection of useful functions to estimate some atmospheric state variables of the marine boundary layer are available in the module **mod_thermo** (mod_thermo.f90).  
Example for computing SSQ of Eq.(1) out of the SST and the SLP:

              PROGRAM TEST_THERMO
                  USE mod_const
                  USE mod_thermo
                  ...

                  SSQ(:,:) = 0.98*q_sat(SST, SLP)
                  ...
              END PROGRAM TEST_THERMO

</div>

<div style="margin: 35px 40px 0px 40px; font-family: Courier; font-size: 24px;width: 960px;" align="left">**> Acknowledging AeroBulk**</div>

<div style="margin: 10px 40px 0px 40px; font-family: Courier; font-size: 20px;width: 960px;" align="justify">_To acknowledge/reference AeroBulk in your scientific work, please cite:_  
Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically significant effects of some approximations in the bulk parameterizations of turbulent air-sea fluxes. _J. Phys. Oceanogr._, 10.1175/JPO-D-16-0169.1. [ [DOI](http://dx.doi.org/10.1175/JPO-D-16-0169.1) ]  
</div>

<center>

<div style="margin: 35px 40px 0px 60px; font-family: Courier New; font-size: 14px;width: 960px;" align="center"><u>Contact</u>  

[**AeroBulk**](http://sourceforge.net/projects/aerobulk/) / L. Brodeau / 2016</div>

</center>
