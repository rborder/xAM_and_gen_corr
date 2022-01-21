#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))      ## for cli
suppressPackageStartupMessages(library(Matrix))        ## for checking feasibility
suppressPackageStartupMessages(library(MASS))          ## drawing MVN random variates
suppressPackageStartupMessages(library(RcppArmadillo)) ## fast lm with std errs
suppressPackageStartupMessages(library(stringr))       ## string wrangling
suppressPackageStartupMessages(library(MatchIt))       ## matching mates
options(digits=3)

###############################################################################
############################## GENOTYPE HANDLING ##############################
###############################################################################

## haploid genotypes -> diploid genotypes
hap_to_dip <- function(G) {
  m2 <- ncol(G)
    G[ ,seq(1,m2,2)] + G[,seq(2,m2,2)]
}

## generate effect sizes, initial genotypes
gen_architect <- function(n, m, rho_beta, s2_y0, s2_z0, theshold_MAF = .1) {
    S = diag(sqrt(c(s2_y0,s2_z0)))
    R_beta <- matrix(c(1,rho_beta,rho_beta,1),2,2)
    C_beta = (S %*% R_beta %*% S)
    beta = mvrnorm(m, c(0,0), C_beta)/sqrt(m)
    beta = cbind(c(rep(0,m/2),rnorm(m/2,0,S[1,1])),
                 c(rnorm(m/2,0,S[2,2]),rep(0,m/2)))/sqrt(m/2)
    q = runif(m,theshold_MAF,1-theshold_MAF)
    q_hap = rep(q,each=2)
    G_raw = sapply(q_hap, function(p) rbinom(n,1,p)) ## n x 2*m
    dip_sd = rep(apply(hap_to_dip(G_raw),2,sd),each=2)
    G_std = sapply(1:(2*m), function(j) (G_raw[,j] - mean(G_raw[,j]))/dip_sd[j])
    return(list(G=G_std, beta=beta, MAF=q))
}

## architecture -> phenotypes, liabilities, noise
compute_pheno <- function(sim) {
    X <- hap_to_dip(sim$data$G)
    sim$data$liability <- X %*% sim$data$beta
    sim$data$noise <- cbind(rnorm(sim$n,0,sqrt(1-sim$s2_y0)),
                   rnorm(sim$n,0,sqrt(1-sim$s2_z0)))
    sim$data$yz <- sim$data$liability + sim$data$noise
    return(sim)
}

## (parent haploid genotypes as nx2m array, mate indices) -> offspring haploid genotypes
## note: halves population size at each generation
## test via:
## matings =rbind(1:13,26:14)
## G = t(sapply(letters,function(x) paste0(x,'_',rep(1:400,each=2),'_',rep(0:1,400))))
meiosis <- function(G, matings, rfreq) {
    n = 2*(nrow(G)%/%2)
    G2 <- G
    m  <- ncol(G) / 2
    ## recombination
    if (length(rfreq)==1) {
        rfreq <- c(.5, rep(rfreq,m-1))
        tmp=seq(1,length(rfreq), length(rfreq)%/%20)
        rfreq[tmp[1:20]] <- .5
    }
    ## meiosis
    G2[ , seq(1,2*m,2)] <- t(apply(G[matings[1,],],1, function(row) row[cumsum(rbinom(m,1,rfreq) )%%2+ seq(1,2*m,2)])) 
    G2[ , seq(2,2*m,2)] <- t(apply(G[matings[2,],],1, function(row) row[cumsum(rbinom(m,1,rfreq) )%%2+ seq(1,2*m,2)])) 
    return(G2)
}


###############################################################################
############################## MATING FUNCTIONS ###############################
###############################################################################
## hyper parameter for matching
matching_multiplier <- function(ryy,rzz,ryz,syz){
    (0.004824 + 0.451851* ryy + 0.369544* ryz + 0.44378* rzz - 
     0.11423* syz)/(0.429722* ryy + 0.355511* ryz + 0.421527* rzz - 
                    0.101386* syz)
}

## match mates
mate_matching_scaled <- function(yz, par) {
    r_yy <- par$r_yy; r_zz <- par$r_zz
    n = 2*(nrow(yz)%/%2)
    ## choose cross mate cross trait corr
    s_yz <- cor(yz[,1],yz[,2])
    if (is.null(par$r_yz)) {
        r_yz <- sqrt(s_yz**2 * abs(r_zz*r_yy))* sign(r_zz*r_yy)
    } else {
        r_yz <- par$r_yz
    }
    mm <- matching_multiplier(r_yy,r_zz,r_yz,s_yz)
    r_yz <- mm*r_yz; r_yy <- mm*r_yy
    r_zz <- mm*r_zz
    target_corr <- matrix(c(1.00,s_yz,r_yy,r_yz,
                            s_yz,1.00,r_yz,r_zz,
                            r_yy,r_yz,1.00,s_yz,
                            r_yz,r_zz,s_yz,1.00), 4, 4)
    ## ensure target correlation structure is feasible
    pdcheck <- Matrix::nearPD(target_corr,corr = T)
    target_corr <- pdcheck$mat 
    ## dummy variates with target distribution and match
    dummy <- apply(mvrnorm(n, rep(0,4), target_corr),2,scale)
    yz <- apply(yz,2,scale)
    yz <- rbind(yz,yz)
    stacked <- as.data.frame(rbind(yz, rbind(dummy[,1:2],dummy[,3:4])))
    names(stacked) <- c('y', 'z')
    stacked$yz <- stacked$y*stacked$z
    stacked$yzdiff <- stacked$y-stacked$z
    stacked$sex <- rep(0:1, each=2*n)
    matches <- matchit(sex~y+z+yz+yzdiff, stacked, #m.order='random',
                       distance = 'mahalanobis',
                       verbose=TRUE)
    m1_inds <- as.integer(matches$match.matrix[1:(n)]) %% n
    m1_inds[m1_inds==0] <- n
    m2_inds <- as.integer(matches$match.matrix[(n+1):(2*n)]) %% n
    m2_inds[m2_inds==0] <- n
    return(rbind(m1_inds,m2_inds))
}

## confirm cross-mate correlation structure is feasible
check_feasibility <- function(ryy,rzz,ryz,syz) {
    corr <- matrix(c(1,syz,ryy,ryz,
                     syz,1,ryz,rzz,
                     ryy,ryz,1,syz,
                     ryz,rzz,syz,1),4,4)
    min(eigen(corr)$values)>=0
}

## helper for testing
check_mating <- function(mfun,yz=test, par=test_par) {
    couplings <- mfun(yz,par)
    mate_phen <- cbind(yz[couplings[1,],],yz[couplings[2,],])
    cor(mate_phen)
}


###############################################################################
########################### SIMULATION / ESTIMATION ###########################
###############################################################################

## estimate betas give phenotypes, stdized genotypes
estimate_prs_weights <- function(y, Z) {
    ## add intercept
    m = ncol(Z)
    if (abs(sd(y)-1) > 1e-5) y <- (y)/sd(y) ## stdize needed
    Zstd <- cbind(rep(1,nrow(Z)), Z)
    tmp <- sapply(1:m, function(j){
        fit <- (fastLmPure(Zstd[ ,c(1,j+1)],y))
        c(fit$coefficients[2], fit$coefficients[2]/fit$stderr[2])})
    return(tmp[1,]) # return point estimates only
}


## print simulation results
sim.print <- function(sim) {
    print(paste(c("############## GENERATION ", str_wrap(sim$generation,2), " ##############"), collapse=''))
    print(""); print(paste(c("Population consists of ", sim$n, " individuals"),collapse=''))
    print(""); print("Phenotypic covariance matrix:")
    print(cov(sim$data$yz))
    print(""); print("Genetic covariance matrix:")
    print(cov(sim$data$liability))
    if (!is.null(sim$data$matings)) {
        print(""); print("Cross-mate phenotypic covariance matrix:")
        print(cov(cbind(sim$data$yz[sim$data$matings[1,],],sim$data$yz[sim$data$matings[2,],])))
    }
    print("");print("")
}

## estimates tr (X %*% t(X) %*% X %*% t(X)) for faster HE
rtrace_K2 <- function(X, l = 500) { 
    u <- replicate(l, rnorm(nrow(X))) ## random probe vectors
    W <- t(u) %*% (X %*% (t(X) %*% (X %*% (t(X) %*% u))))
    sum(diag(W))/l
}

## compute new phenotypes and optionally perform HE regression
sim.update <- function(sim) {
    ii <- sim$generation + 1
    sim$res$generation[ii] <- sim$generation
    ## convert haploid genotypes to standardized diploid genotypes
    G <- hap_to_dip(sim$data$G)
    G_std <- apply(G, 2, scale)
    y_std <- scale(sim$data$yz[,1])
    z_std <- scale(sim$data$yz[,2])
    ## compute true parameters
    ell_y <- G%*%sim$data$beta[,1]
    ell_z <- G%*%sim$data$beta[,2]
    sim$res$h2_y_true[ii] <- var(ell_y)/var(sim$data$yz[,1])
    sim$res$h2_z_true[ii] <- var(ell_z)/var(sim$data$yz[,2])
    sim$res$rho_ell_true[ii]  <- cor(ell_y,ell_z)
    if (sim$HE) { ## perform HE regression
        n = nrow(G)
        m = ncol(G)
        Ky1 <- G_std%*% (t(G_std) %*% y_std)/m
        Ky2 <- G_std%*% (t(G_std) %*% z_std)/m
        trK2 <- rtrace_K2(G_std/sqrt(m))
        denom <- trK2 - n
        sg1 = (t(y_std) %*% Ky1 - t(y_std)%*% y_std)
        sg2 = (t(z_std) %*% Ky2 - t(z_std)%*% z_std)
        covg = (t(y_std) %*% Ky2 - t(y_std)%*% z_std)
        rg = covg / sqrt(sg1*sg2)
        sim$res$h2_y_HE[ii] <- sg1/ denom
        sim$res$h2_z_HE[ii] <- sg2/ denom
        sim$res$rg_HE[ii]   <- rg
    }
    if (sim$PRS) { ## compute GWAS effects + PRS
        beta_hat_y <- estimate_prs_weights(y_std, G_std)
        beta_hat_z <- estimate_prs_weights(z_std, G_std)
        beta_y <- sim$data$beta[,1]
        beta_z <- sim$data$beta[,2]
        sim$res$beta_bias_y[ii] <- mean(sign(beta_y)*(beta_hat_y-beta_y))
        sim$res$beta_bias_z[ii] <- mean(sign(beta_z)*(beta_hat_z-beta_z))
        sim$res$beta_bias_scaled_y[ii] <- mean(abs(beta_hat_y-beta_y))*sqrt(m/sim$s2_y0)
        sim$res$beta_bias_scaled_z[ii] <- mean(abs(beta_hat_z-beta_z))*sqrt(m/sim$s2_z0)
        sim$res$rho_ell_prs[ii] <- cor(G_std%*%beta_hat_y,G_std%*%beta_hat_z)
        sim$res$rho_beta_prs[ii] <- cor(beta_hat_y, beta_hat_z)
    }
   return (sim)
}

## construct and run simulation
simulation <- function(ngen=1, sim = list(n=4000, m=1000, rho_beta=.0, s2_y0=.5,
                                          s2_z0=.5, r_yy=.5, r_zz=0, rfreq=0.5,
                                          data = NULL, PRS=FALSE, HE=FALSE, rmate=0),
                       run_counter = 1, verbose = TRUE, mating_scheme = 'geomean'){
    set.seed(run_counter)
    if (is.null(sim$data)) {
        ## initial genetic architecture
        sim$data = with(sim, gen_architect(n, m, rho_beta, s2_y0, s2_z0))
        sim$data$matings = NULL
        sim$generation = 0
      if (mating_scheme=='xmate') mating_scheme <- sim$r_yz
        ## initialize results
        sim$res <- within(data.frame(generation = rep(NA,ngen+1)), {
            n <- sim$n; m <- sim$m; rho_beta <- sim$rho_beta; s2_y0 <- sim$s2_y0;
            s2_z0 <- sim$s2_z0; r_yy <- sim$r_yy; r_zz <- sim$r_zz;
            h2_y_HE <- h2_z_HE <- rg_HE <- NA
            beta_bias_y <- beta_bias_z <- rho_ell_prs  <- NA
            beta_bias_scaled_y <- beta_bias_scaled_z  <- NA
            rmate <- sim$rmate
            h2_y_true <- h2_z_true <- rho_ell_true <- rho_beta_prs <- NA
            rx_yy <- rx_yz <- rx_zz <- rw_yz <- NA
            rfreq <- sim$rfreq
            seed <- run_counter
            rho_beta_true <- cor(sim$data$beta[,1],sim$data$beta[,2])
            mating_scheme <- mating_scheme
        })
    }
    if (sim$rmate > 0) {
      ## iterate over random mating generations
      for (RGEN in 1:sim$rmate) {
        ## not implemented
      }
    }
    ## iterate over assortative mating generations
    for (i in 1:(ngen+1)) {
        ## compute phenotypes
        sim <- compute_pheno(sim)
        sim <- sim.update(sim)
      tmp_par = list(r_yy = sim$r_yy, r_zz = sim$r_zz)
      ## if using specified xtrait xmate corr
        if (is.numeric(mating_scheme)) {
            tmp_par$r_yz <- sim$r_yz
        }
      ## if using default xtrait xmate corr
      if ( mating_scheme=='geomean' ) { 
            tmp_par <- c(tmp_par, r_yz = sign(sim$r_zz*sim$r_yy)* sqrt(abs(sim$r_zz*sim$r_yy)))
      }
      ## select mates
        sim$data$matings = mate_matching_scaled(yz=sim$data$yz,
                                                par = tmp_par)
        mate_phen <- cbind(sim$data$yz[sim$data$matings[1,],],
                           sim$data$yz[sim$data$matings[2,],])
      ## compute sample xmate correlations
        mate_phen_corr <- cor(mate_phen)
        phen_corr <- cor(sim$data$yz)
        sim$res$rw_yz[i] <- phen_corr[2,1]
        sim$res$rx_yy[i] <- mate_phen_corr[3,1]
        sim$res$rx_zz[i] <- mate_phen_corr[4,2]
        sim$res$rx_yz[i] <- (mate_phen_corr[3,2] + mate_phen_corr[4,1])/2
        sim$n <- ncol(sim$data$matings)
      ## perform meiosis to generate next generation genotypes
        sim$data$G <- meiosis(sim$data$G,sim$data$matings,sim$rfreq)
        if (verbose) sim.print(sim)
        sim$res$generation[i] <- sim$generation
        sim$generation <- sim$generation+1
    }
    return(sim)
}



###############################################################################
##################################### CLI #####################################
###############################################################################
parser <- ArgumentParser()

parser$add_argument("-n", type="integer",
                    help="Founder population size",
                    default=4000)
parser$add_argument("-m", type="integer",
                    help="Number of causal variants",
                    default=1000)
parser$add_argument("-g", type="integer", default=3,
                    help="Number of generations to run")
parser$add_argument("--nsim", type="integer", default=1,
                    help="Number of simulations to run")
parser$add_argument("--rho_beta", type="double", default = 0, 
                    help="Effect size correlation")
parser$add_argument("--rxyy", type="double", default=.5, 
                    help="Corr(Y*, Y**)")
parser$add_argument("--rxzz", type="double", default=.5, 
                    help="Corr(Z*, Z**)")
parser$add_argument("--rxyz", type="double", default=.5, 
                    help="Corr(Y*, Z**)")
parser$add_argument("--h2_y", type="double", default=.5,
                    help="Panmictic h2 of y")
parser$add_argument("--h2_z", type="double", default=.5, 
                    help="Panmictic h2 of z")
parser$add_argument("--mating-scheme", type='character', default = 'xmate',
                    help="How to decide cross mate correlations")
parser$add_argument("--rfreq", type="double", default = 0.5,
                    help="Probability of recombination event between contiguous loci. Controls LD structure")
parser$add_argument("-o", "--output", type="character", default='',
                    help="Output destination")
parser$add_argument("--printfreq", type="integer", default = 5,
                    help="Display progress every ____ iterations")
parser$add_argument("--rmate", type="integer", default = 0,
                    help="number of generations of random mating prior to start")
parser$add_argument("--countfrom", '--seed', type="integer", default=sample.int(5000,1),
                    help="Where to start indexing runs. Serves as random seed")
args <- parser$parse_args()

## choose default output file if not provided
if (args$output=='') {args$output <- paste0(c('test_',args$countfrom,'.out'), collapse='')}

## make sure supplied xmate correlations are feasible
if (!check_feasibility(args$rxyy,args$rxzz,args$rxyz,0)) {
    message('infeasible correlation structure. exiting...')
    quit(save='no')
}
## correct out of bounds values
if (abs(args$h2_y)>1) {
args$h2_y <- sign(args$h2_y)
}
if ((args$h2_y)<0) {
args$h2_y <- 0
}
if (abs(args$h2_z)>1) {
args$h2_z <- sign(args$h2_z)
}
if ((args$h2_z)<0) {
args$h2_z <- 0
}
if (abs(args$rxzz)>1) {
args$rxzz <- sign(args$rxzz)
}
if (abs(args$rxyy)>1) {
args$rxyy <- sign(args$rxyy)
}


## run simulations
run_counter=(0:(args$nsim-1))+args$countfrom
results <-
    lapply(1:args$nsim, function(JJ) {
        message(paste(mapply(function(x,y) paste(x,y, sep=':'), names(args), args),collapse ='  '))
                simulation(ngen=args$g,
                           sim = list(n=args$n, m=args$m, rho_beta=args$rho_beta,
                                      s2_y0=args$h2_y, s2_z0=args$h2_z,
                                      r_yy=args$rxyy, r_zz=args$rxzz,r_yz=args$rxyz,
                                      rfreq=args$rfreq, data = NULL,
                                      HE=TRUE,PRS=TRUE,
                                      rmate=args$rmate),
                           run_counter=run_counter[JJ], verbose = FALSE,
                           mating_scheme=args$mating_scheme)$res
    })
## aggregate results
output <- do.call(rbind.data.frame,results)

## print clean summary
tidy_summary <- aggregate(output[setdiff(names(output),c('generation','mating_scheme'))],
          output[c('generation','mating_scheme')],mean)
cat('\nSummary\n\n')
print(tidy_summary)

## save output
write.table(output, args$output, quote=F, row.names=F, sep=',', col.names=F, append=T)
