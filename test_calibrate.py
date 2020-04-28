import spotpy
from optimizer import spotpy_setup

optimizer = spotpy_setup()

sampler = spotpy.algorithms.mc(optimizer, dbname='MC', dbformat='csv')
# sampler = spotpy.algorithms.mle(spotpy_setup,dbname='MLE',dbformat='csv')
# sampler = spotpy.algorithms.lhs(spotpy_setup,dbname='LHS',dbformat='csv')
# sampler = spotpy.algorithms.sceua(spotpy_setup,dbname='SCEUA',dbformat='csv')
# sampler = spotpy.algorithms.demcz(spotpy_setup,dbname='DE-MCz',dbformat='csv')
# sampler = spotpy.algorithms.sa(spotpy_setup,dbname='SA',dbformat='csv')
# sampler = spotpy.algorithms.rope(spotpy_setup,dbname='ROPE',dbformat='csv')


algorithms=['MC']#,'LHS','MLE','MCMC','SCE-UA','SA','DE-MCz','ROPE']
results=[]
for algorithm in algorithms:
    sampler.sample(10000)
    results.append(sampler.getdata)


