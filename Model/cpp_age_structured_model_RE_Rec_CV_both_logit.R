setwd("D:\\Project\\Thysis\\Model")
library(TMB);

age_structured_model_RE_Rec_CV_both_logit <- '
//State-space age-structured model
//Author: Kim as of Feb 17, 2026
//Random-effects: Rec, Ft; (Fixed CV for BOTH Yield and CPUE)
  
#include <TMB.hpp>
#include <iostream>
  
// pass missing values
template<class Type>
bool isNA(Type x) {
return R_IsNA(asDouble(x));
}
  
// square
template <class Type>
Type square(Type i) {
  return i*i; 
}
  
template <class Type>
Type logDiriMultinom_f(int nages, Type SampleSize, Type theta, vector<Type> prop_obs, vector<Type> prop_pred) {
  Type zero = Type(0.0);
  Type one = Type(1.0);
  
  vector<Type> freq_obs(nages);
  vector<Type> freq_pred(nages); 
  freq_obs = SampleSize * prop_obs;
  freq_pred = SampleSize * prop_pred;
  
  Type ll = zero;
  
  // Dirichlet-Multinomial Log-likelihood
  ll += lgamma(SampleSize + one) + lgamma(SampleSize * theta) - lgamma(SampleSize + SampleSize * theta);
  
  for(int a = 0; a < nages; a++) {
      ll += -lgamma(freq_obs(a) + one) 
            + lgamma(freq_obs(a) + freq_pred(a) * theta + Type(1.0e-10)) 
            - lgamma(freq_pred(a) * theta + Type(1.0e-10));
  }
  
  return ll; 
}
  
//objective function
template<class Type>
Type objective_function<Type>::operator() () {
  
  //Data
  DATA_INTEGER(nyears); // the number of years
  DATA_INTEGER(nages);  // the number of ages
  DATA_MATRIX(WAA_ini); // Weight at age
  DATA_MATRIX(CAA_obs); // age composition
  DATA_VECTOR(yield);
  DATA_VECTOR(CPUE);
  DATA_VECTOR(lambda);  //Weight Parameter
  DATA_SCALAR(Spawn_month); //Spawning month
  DATA_VECTOR(range_q);
  DATA_SCALAR(M);       // input value
  DATA_VECTOR(maturate);
  DATA_VECTOR(Fref);
  DATA_VECTOR(SS_age);
  DATA_VECTOR(musig_logit_h); // 
  DATA_SCALAR(ratio_female);  // ratio_female
  DATA_SCALAR(CV_yield);
  DATA_SCALAR(CV_cpue);
  
  Type q_lower = range_q(0);
  Type q_upper = range_q(1);
  
  // Recruitments (Random-effects) 
  PARAMETER(log_Rec_fir); 
  Type Rec_fir = exp(log_Rec_fir);
  
  // Selectivity 
  PARAMETER(log_a50);                     //gear selectivity parameter;
  PARAMETER(log_a95);                 //gear selectivity parameter;
  Type a50=exp(log_a50);              //logistic selectivity;
  Type a95=exp(log_a95);          //logistic selectivity;
  
  // biomass index 
  PARAMETER(logit_q);  //one scalar because of one survey  

  // Recruits
  PARAMETER(logR0);      // Unexploited Recruitment
  PARAMETER(logit_h);        //Steepness
  
  // Random effect; variance of process error; Rec; Ft; 
  PARAMETER_VECTOR(logRec);  
  PARAMETER(log_sig_logRec); 
  PARAMETER(logF1);
  Type F1 = exp(logF1);
  PARAMETER_VECTOR(logFt); 
  PARAMETER(log_sig_logF); 
  
  // Dirichlet
  PARAMETER(logtheta); 
  Type theta = exp(logtheta);

  //Derived parameters
  vector<Type> Ft=exp(logFt); // length 25
  vector<Type> Sel(nages); 
  Sel.setZero(); 
  vector<Type> M_y(nyears); 
  M_y.setZero(); 
  matrix<Type> F_ta(nyears,nages); 
  F_ta.setZero(); 
  matrix<Type> Z_ta(nyears,nages); 
  Z_ta.setZero(); 
  matrix<Type> ExpZ_ta(nyears,nages); 
  ExpZ_ta.setZero(); 
  
  Type one = Type(1.0);
  
  matrix<Type> WAAcatchall(nyears, nages); 
  matrix<Type> WAAssb(nyears, nages);
  matrix<Type> index_WAA(nyears, nages);
  matrix<Type> NAA(nyears,nages); 
  NAA.setZero(); 
  vector<Type> Biomass_temp(nyears); 
  Biomass_temp.setZero();   
  matrix<Type> CAA_pred(nyears,nages);     
  CAA_pred.setZero(); 
  matrix<Type> Yield_ta(nyears,nages); 
  Yield_ta.setZero(); 
  matrix<Type> B_ta(nyears,nages); 
  B_ta.setZero(); 
  vector<Type> Yield_hat(nyears); 
  Yield_hat.setZero(); 
  matrix<Type> predcpue_ya(nyears, nages); 
  predcpue_ya.setZero();
  vector<Type> predcpue(nyears); 
  predcpue.setZero();
  
  int indnum =  Fref.size();
  vector<Type> recruits_pred(nyears); 
  recruits_pred.setZero();
  matrix<Type> Spawn_NAA(nyears,nages); 
  Spawn_NAA.setZero(); 
  matrix<Type> Spawn_BAA(nyears,nages); 
  Spawn_BAA.setZero(); 
  vector<Type> SSB(nyears); 
  SSB.setZero(); 
  vector<Type> Rec(nyears); 
  Rec.setZero(); 
  vector<Type> YPR(indnum); 
  YPR.setZero(); 
  vector<Type> SPR(indnum); 
  SPR.setZero(); 
  vector<Type> SPR_T(indnum); 
  SPR_T.setZero(); 
  vector<Type> SPR_0(indnum); 
  SPR_0.setZero(); 
  matrix<Type> Zref(indnum,nages); // size: (nind x nages); 
  Zref.setZero(); 
  
  for (int i=0;i<nyears;i++)  { 
    WAAcatchall.row(i)=WAA_ini.row(i); 
    WAAssb.row(i)=WAA_ini.row(i);     
    index_WAA.row(i)=WAA_ini.row(i); 
  }; 
  
  vector<Type> nll(6); // components of the objective function; negative log likelihood 
  
  //get_selectivity (logistic) 
  for (int j=0;j<nages;j++) {
    Type age = Type(j) + 1.0;
    Sel(j)= one/(one+exp(Type(-1.0)*log(19)*((age-a50)/(a95-a50)))); 
  };
  
  //get_mortality_rates; 
  for(int i=0; i<nyears; i++){ 
    M_y(i)=M; 
  }; 
  for (int i=0; i<nyears; i++){ 
    for (int j=0; j<nages; j++){ 
      F_ta(i,j) = Sel(j)*Ft(i); 
      Z_ta(i,j) = M_y(i)+F_ta(i,j); 
      ExpZ_ta(i,j) = exp(Type(-1.0)*Z_ta(i,j)); 
    }; 
  }; 
  
  Type R0 = exp(logR0);
  Type h = Type(0.2001) + Type(0.7999) / (Type(1.0) + exp(-logit_h));
  Type SR_alpha;
  Type SR_beta;
  Type SSB0;
  
  // SPR0
  vector<Type> SPR_0_vec(nyears);
  Type M_val = M_y(nyears - 1); 
  
  for (int y = 0; y < nyears; y++) {
    SPR_0_vec(y) = 
        (WAAssb(y, 0) * maturate(0) * ratio_female * exp(-(Spawn_month / 12.0) * M_val) + 
        exp(-M_val) * WAAssb(y, 1) * maturate(1) * ratio_female * exp(-(Spawn_month / 12.0) * M_val) + 
        exp(-M_val * 2.0) * WAAssb(y, 2) * maturate(2) * ratio_female * exp(-(Spawn_month / 12.0) * M_val) + 
        exp(-M_val * 3.0) * WAAssb(y, 3) * maturate(3) * ratio_female * exp(-(Spawn_month / 12.0) * M_val) + 
        exp(-M_val * 4.0) * WAAssb(y, 4) * maturate(4) * ratio_female * exp(-(Spawn_month / 12.0) * M_val) + 
        (exp(-M_val * 5.0)) / (Type(1.0) - exp(-M_val)) * WAAssb(y, 5) * maturate(5) * ratio_female * exp(-(Spawn_month / 12.0) * M_val))/Type(1000);
  }

  int y0 = 0;
  // Steady-state
  NAA(y0, 0) = exp(log_Rec_fir); 
  for (int a = 1; a < nages - 1; a++) {
    NAA(y0, a) = NAA(y0, a - 1) * ExpZ_ta(y0, a - 1);
  }
  NAA(y0, nages - 1) = NAA(y0, nages - 2) * ExpZ_ta(y0, nages - 2) / (Type(1.0) - ExpZ_ta(y0, nages - 1));
  
  for (int a = 0; a < nages; a++) {
      Spawn_NAA(y0, a) = NAA(y0, a) * maturate(a) * ratio_female * exp(-(Spawn_month / 12.0) * Z_ta(y0, a));
      Spawn_BAA(y0, a) = Spawn_NAA(y0, a) * WAAssb(y0, a) / Type(1000.0);
  }
  SSB(y0) = Spawn_BAA.row(y0).sum();
  
  // first recruits
  recruits_pred(y0) = (Type(4.0) * h * SSB(y0) * R0) / (R0 * SPR_0_vec(y0) * (Type(1.0) - h) + SSB(y0) * (Type(5.0) * h - Type(1.0)));
  
  // Cohort loop (Year 2 to nyears)
  for (int y = 1; y < nyears; y++) {
      NAA(y, 0) = exp(logRec(y));
      for (int a = 1; a < nages; a++) {
          if (a < nages - 1) {
               NAA(y, a) = NAA(y - 1, a - 1) * ExpZ_ta(y - 1, a - 1);
          } else {
               NAA(y, a) = NAA(y - 1, a - 1) * ExpZ_ta(y - 1, a - 1) + NAA(y - 1, a) * ExpZ_ta(y - 1, a);
          }
      }

      for (int a = 0; a < nages; a++) {
          Spawn_NAA(y, a) = NAA(y, a) * maturate(a) * ratio_female * exp(-(Spawn_month / 12.0) * Z_ta(y, a));
          Spawn_BAA(y, a) = Spawn_NAA(y, a) * WAAssb(y, a) / Type(1000.0);
      }
      SSB(y) = Spawn_BAA.row(y).sum();

      recruits_pred(y) = (Type(4.0) * h * SSB(y) * R0) / (R0 * SPR_0_vec(y) * (Type(1.0) - h) + SSB(y) * (Type(5.0) * h - Type(1.0)));
  }
  Rec = NAA.col(0);
  Type aver_SPR_0 = SPR_0_vec.sum() / nyears;
  SR_alpha = Type(4.0) * h * R0 / ((Type(5.0) * h) - Type(1.0));
  SR_beta  = R0 * aver_SPR_0 * (Type(1.0) - h) / (Type(5.0) * h - Type(1.0));
  SSB0     = aver_SPR_0 * R0;
  
  Type q = q_lower + (q_upper - q_lower) / (Type(1.0) + exp(-logit_q));
  vector<Type> Qa(nages);
  Qa.setZero();
  for(int a = 0; a < nages; a++) {
      Qa(a) = q * Sel(a); 
  }

  for (int y = 0; y < nyears; y++) {
      for (int a = 0; a < nages; a++) {
          CAA_pred(y, a) = NAA(y, a) * (F_ta(y, a) / Z_ta(y, a)) * (Type(1.0) - ExpZ_ta(y, a));
          Yield_ta(y, a) = CAA_pred(y, a) * WAAcatchall(y, a);
          Type B_start = NAA(y, a) * WAAcatchall(y, a); 
          B_ta(y, a) = B_start;
          Type BDY = Type(0.5) * (B_start + (B_start - Yield_ta(y, a))); 
          predcpue_ya(y, a) = BDY * Qa(a);
      }
      Yield_hat(y) = Yield_ta.row(y).sum();
      Biomass_temp(y) = B_ta.row(y).sum();
      predcpue(y) = predcpue_ya.row(y).sum() / Type(1000.0);
   }
   
  //Objective function 
  nll.setZero(); 
  
  // Dirichlet multinomial for age-composition 
  for(int y = 0; y < nyears; y++) {
      vector<Type> prop_obs = CAA_obs.row(y) / CAA_obs.row(y).sum();
      vector<Type> prop_pred = CAA_pred.row(y) / CAA_pred.row(y).sum();
      nll(0) -= logDiriMultinom_f(nages, SS_age(y), theta, prop_obs, prop_pred); 
  } 
  nll(0) = lambda(0) * nll(0);
  std::cout << "nll(0): " << nll(0) << std::endl;
  
  // observation error - cv_yield
  Type sig2_yield = log(square(CV_yield) + one);  
  Type sig_yield = sqrt(sig2_yield);
  for(int y = 0; y < nyears; y++) {
      nll(1) -= dnorm(log(yield(y)), log(Yield_hat(y) / Type(1000.0)), sig_yield, true);
  }
  nll(1) = lambda(1) * nll(1);
  std::cout << "nll(1): " << nll(1) << std::endl;
  
  // observation error - cv_cpue
  Type sig2_cpue = log(square(CV_cpue) + one);  
  Type sig_cpue = sqrt(sig2_cpue);
  for(int y = 0; y < nyears; y++) {
      nll(2) -= dnorm(log(CPUE(y)), log(predcpue(y)), sig_cpue, true);
  }
  nll(2) = lambda(2) * nll(2);
  std::cout << "nll(2): " << nll(2) << std::endl;
  
  // fishing mortality; penalized 
  Type sig_logF = exp(log_sig_logF);
  nll(3) -= dnorm(logFt(0), logF1, sig_logF, true);  
  for(int y = 1; y < nyears; y++) {
      nll(3) -= dnorm(logFt(y), logFt(y - 1), sig_logF, true);
  } 
  nll(3) = lambda(3) * nll(3);
  std::cout << "nll(3): " << nll(3) << std::endl;
  
  //process of recruits as random effects
  Type sig_logRec = exp(log_sig_logRec);
  nll(4) -= dnorm(logRec(0), log(NAA(0,0)), sig_logRec, true);
  for(int y = 1; y < nyears; y++) {
      nll(4) -= dnorm(logRec(y), log(recruits_pred(y - 1)), sig_logRec, true);
  } 
  nll(4) = lambda(4) * nll(4);
  std::cout << "nll(4): " << nll(4) << std::endl;
  
  //Steepness
  Type h_prior = musig_logit_h(0);  
  Type h_logit_sd = musig_logit_h(1); 
  Type mu_logit = log(((h_prior - Type(0.2001)) / Type(0.7999)) / (Type(1.0) - (h_prior - Type(0.2001)) / Type(0.7999)));
  nll(5) -= dnorm(logit_h, mu_logit, h_logit_sd, true);
  nll(5) = lambda(5) * nll(5);
  std::cout << "nll(5): " << nll(5) << std::endl;
  
  Type jnll=nll.sum();
  
  /////////////////////////////////////// YPR & SPR ///////////////////////////////////////// 

  for(int ind = 0; ind < indnum; ind++) { 
    for(int j = 0; j < nages; j++) { 
      Zref(ind, j) = Fref(ind) * Sel(j) + M_y(nyears - 1); 
    } 
  }

  // YPR
  for(int ind = 0; ind < indnum; ind++) { 
    YPR(ind) = (WAAcatchall(nyears-1,0)*(Fref(ind)*Sel(0)/Zref(ind,0))*(Type(1.0)-exp(-Zref(ind,0))) + 
               exp(-Zref(ind,0))*WAAcatchall(nyears-1,1)*(Fref(ind)*Sel(1)/Zref(ind,1))*(Type(1.0)-exp(-Zref(ind,1))) + 
               exp(-Zref(ind,0)-Zref(ind,1))*WAAcatchall(nyears-1,2)*(Fref(ind)*Sel(2)/Zref(ind,2))*(Type(1.0)-exp(-Zref(ind,2))) + 
               exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2))*WAAcatchall(nyears-1,3)*(Fref(ind)*Sel(3)/Zref(ind,3))*(Type(1.0)-exp(-Zref(ind,3))) + 
               exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2)-Zref(ind,3))*WAAcatchall(nyears-1,4)*(Fref(ind)*Sel(4)/Zref(ind,4))*(Type(1.0)-exp(-Zref(ind,4))) + 
               (exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2)-Zref(ind,3)-Zref(ind,4)))/(Type(1.0)-exp(-Zref(ind,5)))*WAAcatchall(nyears-1,5)*(Fref(ind)*Sel(5)/Zref(ind,5))*(Type(1.0)-exp(-Zref(ind,5))))/Type(1000); 
  } 

  // SPR_T 
  for(int ind = 0; ind < indnum; ind++) { 
    SPR_T(ind) = (WAAssb(nyears-1,0)*maturate(0)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,0)) + 
                 exp(-Zref(ind,0))*WAAssb(nyears-1,1)*maturate(1)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,1)) + 
                 exp(-Zref(ind,0)-Zref(ind,1))*WAAssb(nyears-1,2)*maturate(2)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,2)) + 
                 exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2))*WAAssb(nyears-1,3)*maturate(3)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,3)) + 
                 exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2)-Zref(ind,3))*WAAssb(nyears-1,4)*maturate(4)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,4)) + 
                 (exp(-Zref(ind,0)-Zref(ind,1)-Zref(ind,2)-Zref(ind,3)-Zref(ind,4)))/(Type(1.0)-exp(-Zref(ind,5)))*WAAssb(nyears-1,5)*maturate(5)*ratio_female*exp(-(Spawn_month/12.0)*Zref(ind,5)))/Type(1000); 
  }

  for(int ind = 0; ind < indnum; ind++) { 
    SPR_0(ind) = SPR_0_vec(nyears - 1); 
    SPR(ind) = SPR_T(ind) / SPR_0(ind); 
  }
  
  vector<Type> SSB_F(indnum);
  SSB_F.setZero();
  vector<Type> Rec_F(indnum);
  Rec_F.setZero();
  vector<Type> Yield_F(indnum);
  Yield_F.setZero();
 
  for(int ind = 0; ind < indnum; ind++){
    Type ssb_temp = SR_alpha * SPR_T(ind) - SR_beta;
    if(ssb_temp < Type(0.0)) {
      ssb_temp = Type(0.0);
    }
    SSB_F(ind) = ssb_temp;
    
    if(SPR_T(ind) > Type(0.0) && ssb_temp > Type(0.0)) {
      Rec_F(ind) = ssb_temp / SPR_T(ind);
    } else {
      Rec_F(ind) = Type(0.0);
    }
    Yield_F(ind) = YPR(ind) * Rec_F(ind);
  }

  Type MSY = max(Yield_F);
  
  // Report
  REPORT(nll); 
  REPORT(jnll); 

  REPORT(yield);            
  REPORT(Yield_hat);        

  REPORT(CPUE);            
  REPORT(predcpue);        

  REPORT(CAA_obs);         
  REPORT(CAA_pred);        
  
  REPORT(NAA);            
  REPORT(B_ta);          
  REPORT(Biomass_temp);   

  REPORT(Sel);             
  REPORT(F_ta);            
  REPORT(q);             
  
  REPORT(h);
  REPORT(sig_logRec);
  
  REPORT(sig_yield);
  REPORT(sig_cpue);
  
  ADREPORT(Rec_fir);
  ADREPORT(a50);
  ADREPORT(a95);
  ADREPORT(q);
  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(sig_logRec);
  ADREPORT(F1);
  ADREPORT(sig_logF);
  ADREPORT(theta);
  
  REPORT(SSB0);            
  REPORT(SR_alpha);        
  REPORT(SR_beta);
  
  REPORT(YPR);             
  REPORT(SPR);             
  REPORT(SPR_T);           
  REPORT(Yield_F);         
  REPORT(MSY);             
  REPORT(Rec_F);           
  REPORT(SSB_F);
  
  return jnll;  
}' 

# 
write(age_structured_model_RE_Rec_CV_both_logit, file="age_structured_model_RE_Rec_CV_both_logit.cpp");  
compile("age_structured_model_RE_Rec_CV_both_logit.cpp"); 
dyn.load(dynlib("age_structured_model_RE_Rec_CV_both_logit"));
