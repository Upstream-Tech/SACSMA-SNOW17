1) Compile the Snow17 module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/snow19
>> f2py -c -m exsnow19.f *.f

2) Compile the SAC-SMA module.
>> cd /Users/grey/workspace/SACSMA-SNOW17/sacsma_source/sac
>> 

3) Compile the old SAC-SMA module
>> f2py -c sacsma.pyf sacsma.f 
