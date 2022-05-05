## NO constraints (all 24 parameters)
full_pars <- list()

## least-constrained model (12 params: 24 → collapse 4 sets of 4 to 1 each

pcsc_pars <- list(c(7, 10, 20, 23),  ## all pc_gain rates
             c(4, 11, 17, 24),  ## all sc_gain rates
             c(2, 5, 15, 18),   ## all pc_loss rates
             c(1, 8, 14, 21))   ## all sc_loss rates

## add constraints: ag gain/loss depends only on pc
pc_pars <- c(pcsc_pars,
          list(c(13, 16), ## ag_gain with pc==0 (sc==0 or 1)
               c(19, 22), ## ag_gain with pc==1 (sc==0 or 1)
               c(3,   6), ## ag_loss with pc==0
               c(9,  12)) ## ag_gain with pc==1
          )

## add constraints: ag gain/loss depends only on sc
sc_pars, <- c(pcsc_pars,
          list(c(13, 19), ## ag_gain with sc==0
               c(16, 22), ## ag_gain with sc==1
               c(3,   9), ## ag_loss with sc==0
               c(6,  12)) ## ag_gain with sc==1
          )
