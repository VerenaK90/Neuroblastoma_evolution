import math

from pyabc.external import R

tumor_id=''

rname = "~/Data/" + tumor_id + "/Neutral_fit_pre_clonal_and_clonal.R"
r = R(rname)

model = r.model("myModel")
distance = r.distance("mySummaryDistance")
observation = r.observation("mySumStatData")


from pyabc import Distribution, RV, ABCSMC 


prior = Distribution(n_clonal=RV("uniform", 0, 10000), mu=RV("uniform", 0.1, 19.9), delta=RV("uniform", 0, 0.99))

from pyabc.sampler import MulticoreEvalParallelSampler

multi_sampler = MulticoreEvalParallelSampler(n_procs=16)

abc = ABCSMC(model, prior, distance, sampler = multi_sampler, population_size = 1000)

import os
from tempfile import gettempdir

db = "Data/" + tumor_id + "/" + tumor_id + ".db"

abc.new(db, r.observation("mySumStatData"))

history = abc.run(minimum_epsilon=0.05, max_nr_populations=25)


