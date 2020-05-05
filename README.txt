This repository is for compiling legacy SAC-SMA fortran code with a Python linker (f2py3). There are three parts to the full build: sacsma, snow17, and duamel (unit hydrograph) - these can be compiled using teh following instructions. The 'sacsma_utilities' function shows how to use these builds. The 'test_sacsma' function allows you to run the builds using NCAR CAMELS data.

1) Compile the Snow17 module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/snow19
>> f2py -c -m exsnow exsnow19.f PACK19.f   SNDEPTH.f  SNEW.f     SNOWPACK.f SNOWT.f    adjc19.f   aeco19.f   aesc19.f  melt19.f   rout19.f   updt19.f   zero19.f

2) Compile the SAC-SMA module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/sac
>> f2py -c -m exsac ex_sac1.f sac1.f

3) Compile the unit hydrgraph router
>> f2py -c sacsma.pyf sacsma.f
~                                                                                                                             
~                              
