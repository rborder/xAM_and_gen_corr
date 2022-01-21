## Code accompanying _Cross-trait assortative mating is widespread and inflates genetic correlation estimates_

### Dependendencies

This code is intended to be run using `R` version 4.0 or newer. The following packages are required:

 - `argparse`
 - `MASS`
 - `MatchIt`
 - `Matrix`
 - `RcppArmadillo`
 - `stringr`

### Usage

Use `Rscript xAM_sim.r --help` to display command line arguments:

```r
usage: xAM_sim.r [-h] [-n N] [-m M] [-g G] [--nsim NSIM] [--rho_beta RHO_BETA]
                 [--rxyy RXYY] [--rxzz RXZZ] [--rxyz RXYZ] [--h2_y H2_Y]
                 [--h2_z H2_Z] [--mating-scheme MATING_SCHEME] [--rfreq RFREQ]
                 [-o OUTPUT] [--printfreq PRINTFREQ] [--rmate RMATE]
                 [--countfrom COUNTFROM]

optional arguments:
  -h, --help            show this help message and exit
  -n N                  Founder population size
  -m M                  Number of causal variants
  -g G                  Number of generations to run
  --nsim NSIM           Number of simulations to run
  --rho_beta RHO_BETA   Effect size correlation
  --rxyy RXYY           Corr(Y*, Y**)
  --rxzz RXZZ           Corr(Z*, Z**)
  --rxyz RXYZ           Corr(Y*, Z**)
  --h2_y H2_Y           Panmictic h2 of y
  --h2_z H2_Z           Panmictic h2 of z
  --mating-scheme MATING_SCHEME
                        How to decide cross mate correlations
  --rfreq RFREQ         Probability of recombination event between contiguous
                        loci. Controls LD structure
  -o OUTPUT, --output OUTPUT
                        Output destination
  --printfreq PRINTFREQ
                        Display progress every ____ iterations
  --rmate RMATE         number of generations of random mating prior to start
  --countfrom COUNTFROM, --seed COUNTFROM
                        Where to start indexing runs. Serves as random seed
```

### Example run

```bash
Rscript xAM_sim.r \
    -n 4000 -m 1000 -g 5 --nsim 1 \
    --rho_beta 0 --h2_y .5 --h2_z .5 \
    --rxyy .5 --rxzz .5 --rxyz .25 \
    --seed 123
```

