## Synthesis of Advice on Priors

Jess' advice was since retracted after further explaination. Sigal's advice fits into John's so basically John's advice should cover it and is as follows:

for PC: loss 10X faster than gain - gain 5X faster than loss

for SC: loss 5X faster than gain - gain 10X faster than loss

for AG: loss 100-1000X faster than gain - gain 10X faster than loss

baseline rates

scale edge lengths so that (sum edge length) == 1: **whole-tree** scaling
  then e.g. a rate of 5 would mean "we expect 5 transitions in this trait over
  the whole tree"
  
scale edge lengths so that (sum edge length) = nspecies; **per-species** scaling
then a rate of 5 would mean "we expect the lineage leading to each species to have
experience 5 transitions"

numbers like 10 transitions per species seem very high
numbers like 10 transitions per *tree* seem very low

we could scale the sum of the branch lengths to 1

upper bound is 10*nspecies
lower bound is 10

and say this range is +/- 3 sigma

see how this compares to the default range in corHMM

