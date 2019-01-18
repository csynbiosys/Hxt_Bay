// RStan model containing the first inducer exchange model for the Toggle Switch contained in the Lugagne paper. Bayesian inference is perfrmed using multiple sets of experimental data instead of just one single set

functions{
  // Function containing the ODE to be used for the inference
  real[] Hxt1(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    // Inputs
    real u_Glu = x_r[1];
    // Parameters
    real v = p[1];
    real K = p[2];
    real h = p[3];
    real g = p[4];
    
    
    //Equations
    real dInd_dt[1];
    dInd_dt[1] = v*(x_r[1])^h/(K^h+x_r[1]^h)-g*y[1];
  //RESULTS
    return dInd_dt;
  }
  
  // Function type vector containing the equations where the root needs to be calculated for the steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    vector[1] alpha;
    // Parameters
    real v = p[1];
    real K = p[2];
    real h = p[3];
    real g = p[4];
    // Equations
    alpha[1] = v/g*(x_r[1]^h/(K^h+x_r[1]^h));
    // Results
    return alpha;
  }
  
}

data {
    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    int stsl[1,m]; // Number of elements at each time series for each series m
    int sts[stslm,m]; // Sampling times for each series m
    real GFPmean[stslm,m]; // estimated observables for tetR+GFP at each sampling time
    real GFPstd[stslm,m]; // standard error for tetR+GFP at each sampling time
    
    // Inputs
    int elm; // Number of rows in the matrices for Glu, half for the inputs and -1 for the total number of events
    int tml; // Maximum length of the rows for the sampling times matrix
    int Nsp[1,m]; // Number of event switching points (including final time) times for each series m
    real ts[tml, m]; // Time series for each serie m
    int tsl[1,m]; // Length of sampling time series per serie m
    real preGlu[1,m]; // Values of inputs for each serie m for the ON incubation 
    real Glu[elm,m]; // Values of inputs at each event for each serie m
    real inputs[(elm),m]; // Input values for each Glu event 
    int evnT[(elm+1),m]; // Event change time points for each series m
    
    // Over night incubation times
    int tonil;
    real toni[tonil];

}

transformed data {
  int nParms = 4; // Number of parameters of the model
  int Neq = 1; // Total number of equations of the model
  int x_i[0]; // Empty x_i object (needs to be defined)
  real x_r[(elm),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq,m]; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP
  real pre[1,m]; // Input values during the ON incubation ordered as IPTG, aTc
  
  for(i in 1:m){
    ivss[1,i] = GFPmean[1,i];
    pre[1,i] = preGlu[1,i];
  };
}

parameters {
    // Parameters to be infered in the model
    real v_raw;
    real K_raw;
    real h_raw;
    real g_raw;
}

transformed parameters {
  // Introduction of the paramemters in an indexed object
  real theta[nParms];
  theta[1] = v_raw;
  theta[2] = K_raw;
  theta[3] = h_raw;
  theta[4] = g_raw;
}

model {
  
  // Intermediate parameters
  int i; // Increasing index for the inputs
  vector[1] ing; // Vector that will include the solutio of the algebraic solution for the steady state of the model
  real ssv[tonil,Neq]; // Real that will include the solution of the ODE for the ON incubation (24h)
  real Y0[Neq,m]; // Initial values for the ODEs variables at the first event
  real yhat2[stslm,m]; // Reals that will include the values of RFP and GFP over time separately

  
  // Priors definition (test)
  v_raw ~ lognormal( 9.9379,3.4068);
  K_raw ~ lognormal(1.0161,0.3483);
  h_raw ~ lognormal(7.1019,2.4346);
  g_raw ~ lognormal(0.2235,0.0766);

  
  // Likelihood
  for (j in 1:m){
    real ivst[Neq]; // Initial value of the states 
    real y_hat[(tsl[1,j]),Neq];
    // Calculation of initial guesses
    ing = SteadyState(to_vector(ivss[1,j]), to_vector(theta), pre[1,j], x_i); // Calculation of initial guesses for steady state
    Y0[1,j] = preGlu[1,j];
    Y0[2,j] = ing[1];
    ssv = integrate_ode_bdf(Hxt1, Y0[,j],0,toni,theta,pre[1,j], x_i, 1e-9, 1e-9, 1e7); // ON incubation calculation for the steady state
    
    Y0[,j] = ssv[tonil];
    i = 1;
    
    // Loop (over the number of events) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    for (q in 1:Nsp[1,j]-1){
      
      int itp = evnT[q,j];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop
      // Calculation of the solution for the ODEs where for events that are not the firt one the time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){
        ivst = Y0[,j];
        part1 = integrate_ode_bdf(Hxt1,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i,j]), x_i, 1e-9, 1e-9, 1e7);
      }
      else{
        part1 = integrate_ode_bdf(Hxt1, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i,j]), x_i, 1e-9, 1e-9, 1e7);
      }

      // Modification of the initial state values for the next event
      ivst = part1[lts];
      // Increase index for inputs
      i=i+1;
      
      // Introduction of the result of part1 into the object y_hat
      for (y in (itp+1):(itp+lts)){
        y_hat[(y),]=(part1)[(y-itp),];
      };
    };

    // Likelihood at each sampling time
    for (t in 1:stsl[1,j]){
      yhat2[t,j] = y_hat[(sts[t,j]+1),2];
      GFPmean[t,j] ~ normal(yhat3[t,j],GFPstd[t,j]);
    }

  };

}

