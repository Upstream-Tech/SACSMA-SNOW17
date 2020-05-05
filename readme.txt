1) Compile the Snow17 module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/snow19
>> f2py -c -m exsnow exsnow19.f PACK19.f   SNDEPTH.f  SNEW.f     SNOWPACK.f SNOWT.f    adjc19.f   aeco19.f   aesc19.f  melt19.f   rout19.f   updt19.f   zero19.f

2) Compile the SAC-SMA module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/sac
>> f2py -c -m exsac ex_sac1.f sac1.f

3) Compile the unit hydrgraph router
>> f2py -c sacsma.pyf sacsma.f 
