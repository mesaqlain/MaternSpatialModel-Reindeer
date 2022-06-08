## Written by Murshid Saqlain
## Script for parameter estimation of the ST region in the Reindeer data
## (Skarin and Alam 2017)

# Parameters for creating the mesh
lattice_spacing <- 50 #How fine should the mesh be? 
extension_points <- 8 #How far beyond the observation area should the mesh extend to, to avoid boundary effects

## We use the same mesh as in Regular Lattice case.
mesh <- create.regular.mesh(x.coords, y.coords, lattice_spacing, extension_points)

## Check how the mesh and observations look like
plot_mesh_and_obs(x.coords, y.coords, mesh)

## Matern Parameters for simulating the gaussian responses
range <- c(350, 400)
alpha <- 2

family <- binomial(link=logit) #Response family

### CREATE MODEL FRAME / SETUP MATRICES ###
# Create Model Frame
model_frame <- model.frame(I(Hog>0)~ sqrt(distvindkraft/100) * factor(year_descrip) +
                             factor(smd50upveg3),
                           data=ST_data)

# Design Matrices X
X <- model.matrix(model_frame, data=ST_data)

# Create vector for response, y 
y.response <- as.numeric(model.response(model_frame))
N <- length(y.response) #No. of observations

b_values <- c()
s0_values <- c()
logl_values <- c()

for (j in 1:length(range))
{
  
  # Create Z matrix
  Z <- create.Z(x.coords, y.coords, mesh, range[j], alpha, uid)
  
  #Create the T matrix. (Lee, Neldar, Pawitan 2006, p154)
  T.mat <-rbind(cbind(X,Z),cbind(Matrix(0,nrow=ncol(Z),ncol=ncol(X)),diag(ncol(Z))))
  
  ######## Initial Values ######## 
  #b0  values are taken from Skarin and Alam 2017, 
  #Refer to lines 75-84 for order of the parameters, precip is not included
  b0 <- c(-3.01, 0.48, -0.28, -0.74, 0.08, 0.85, -0.22, -0.27, -0.60) 
  length_b0 <- length(b0)
  phi0 <- 1
  s0 <- 1
  r0 <- c(rep(0,ncol(Z)))
  length_r0 <- length(r0)
  my.fam <- family
  
  Q <- create.Q_matrix_augmented(create.Q_matrix(alpha, range[j], mesh, lattice_spacing), uid)
  #print(proc.time()-start_time)
  
  #message(sprintf("Calculaing determinant of Q"))
  #start_time <- proc.time()
  det.Q <- (determinant(Q,logarithm=TRUE)) ## Decrease computation time?? 
  #print(proc.time()-start_time)
  
  #message(sprintf("Calculating eta."))
  eta <-(X%*%b0) + (Z%*%r0)
  eta <- eta[,1]
  #dim(eta) #2237x1 as expected
  
  HL.correction <- 0
  mu <- family$linkinv(eta) #g(mu) = eta, so mu = g.inv(eta)
  
  colnames(X)
  #[1] "(Intercept)"                         -0.51                                             
  #[2] "sqrt(distvindkraft/100)"              0.01                             
  #[3] "factor(year_descrip)construction"     1.25                     
  #[4] "factor(year_descrip)operation"       -2.21                   
  #[5] "factor(smd50upveg3)Forest"           -1.34                   
  #[6] "factor(smd50upveg3)Clear"             0.09                   
  #[7] "factor(smd50upveg3)Young"            -1.35                   
  #[8] "precip"                              -5.69                       
  #[9] "sqrt(distvindkraft/100):factor(year_descrip)construction"  -0.26
  #[10] "sqrt(distvindkraft/100):factor(year_descrip)operation"     0.14
  
  for(i in 1:200){
    #message(sprintf("Iteration number: %i", i))
    
    #message(sprintf("Calculating working response."))
    #Working Reponse, pg 44 in text Lee, Neldar, Pawitan 2006
    z <- eta + (y.response - mu) / family$mu.eta(eta) 
    z <- z - HL.correction
    #dim(z) #2237x1
    # Sigma_a
    
    #message(sprintf("Calculating weights."))
    wt <- (1/phi0)*as.numeric((family$mu.eta(eta))^2/family$variance(mu))
    
    #message(sprintf("Constructing Sigma_a matrix."))
    Si <- bdiag(diag(wt),Q/s0)
    #dim(Si) # (5273 x 5273) 
    #Makes sense, 2237x2237 for the sigma on top left block (covariance matrix for the observations)
    #and 3036x3036 for D in the bottom right block (covariance matrix for the random effects)
    
    #Construct T (pg 154 in text Lee, Neldar, Pawitan 2006)
    #dim(T) 5273 x 3046 , as we expected. 
    #Because on the left side we have a vector of y and psi(m) that is 5273x1.
    #y is 2237x1 (observations) and psi(m) is 3036x1 (0s for random effects)
    #On the right side, we have T%*%delta. delta is a vector of beta and random effects.
    #beta is 10x1 and random effects are 3036x1, so delta is 3046x1.
    #T should therefore be #5273 x 3046 , as we expected.
    
    #message(sprintf("Taking Cholesky of Sigma_a."))
    #start_time <- proc.time()
    CSI <- chol(Si) #Take Cholesky of Sigma Inverse 
    #print(proc.time()-start_time)
    #dim(CSI) #5273x5273
    
    #message(sprintf("Taking Cholesky of T."))
    T1 <- CSI%*%T.mat #T* 
    #dim(T) # 5273x3046
    
    #message(sprintf("Constructing vector for responses, y_a."))
    ys <- CSI%*%c(as.numeric(z),rep(0,ncol(Z)))
    
    #message(sprintf("Performing QR decomposition of the cholesky of T."))
    QR1 <- qr(T1) #QR decomposition
    
    #message(sprintf("Finding solutions from the QR decomposition."))
    b1 <- qr.coef(QR1,ys)
    
    #message(sprintf("Creating vector for fixed effects."))
    fe <- b1[1:ncol(X)] #Fixed effects
    
    #message(sprintf("Creating vector for random effects."))
    re <- b1[-(1:ncol(X))] #Random effects
    
    #message(sprintf("Updating eta."))
    eta1 <- (X%*%fe) + (Z%*%re)
    
    #message(sprintf("Checking for convergence."))
    if(sum((eta-eta1)^2)<1e-6*sum(eta1^2)) break #Convergence criterion
    
    #s0_opt <-optimize(logL, interval=c(0.00001,100),maximum=TRUE)
    #message(sprintf("Optimizing log likelihood."))
    opt_values <-nlminb(start = c(s0), #SHOULD I CHANGE THIS TO BE AN ARGUMENT??
                        logL_binom, 
                        lower=0.001,
                        upper=10,
                        control = list(rel.tol=1e-6),
                        phi = phi0,
                        re = re,
                        det.Q = det.Q,
                        Q = Q,
                        T.mat = T.mat,
                        mu = mu,
                        eta = eta,
                        family = family,
                        y.response = y.response)
    
    #message(sprintf("Optimized."))
    #s0 <- s0_opt$maximum
    #s0_obj <- s0_opt$objective
    #message(sprintf("Updating marginal variance, fixed effects, and random effects."))
    s0 <- opt_values$par[1]
    #phi0 <- opt_values$par[2]
    b0 <- fe
    r0 <- re
    
    logL_value <- opt_values$objective
    
    eta <- eta1
    eta <- eta[,1]
    
    mu <- family$linkinv(eta) 
    HL.correction <- HL11(fv = mu, Z = Z, family = family, phi = phi0, s = s0)
    
    
    #message(sprintf("%i\t%.4f\t%.4f", i, b0, s0))
    print(b0)
    print(s0)
    print(logL_value)
  }
  b_values <- c(b_values, b0)
  s0_values <- c(s0_values, s0)
  logl_values <- c(logl_values, logL_value)
}

data <- data.frame(b_values, rep(s0_values,each=9), rep(logl_values,each=9))
write.csv(data, 
          "ST_results2.csv", row.names = F, quote = F)


length_b0 <- length(b0)
max_logL_index <- which.max(logl_values)

max_logL <- logl_values[max_logL_index]
s0_best <- s0_values[max_logL_index]
b0_best <-  b_values[(((max_logL_index-1)*length_b0) + 1) : (((max_logL_index-1)*length_b0)+length_b0)]
range_best <- range[max_logL_index]

data_2 <- data.frame(b0_best, rep(s0_best,each=9), rep(max_logL,each=9), rep(range_best,each=9))
write.csv(data_2, 
          "ST_results_best2.csv", row.names = F, quote = F)