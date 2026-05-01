DATA_SECTION

 //DATA
 init_int Nages // number of ages
 init_int Nsizes // number of sizes
 init_int Nyrs // number of years
 init_ivector sizebins(1,Nsizes) // sizebins
 init_ivector ages(1,Nages) // ages
 init_3darray nijk(1,Nyrs,1,Nages,1,Nsizes) // age length sample by year
 init_matrix yjk(1,Nyrs,1,Nsizes) // catch at length by year
 init_matrix xik(1,Nyrs,1,Nages) // age only samples by year
 init_ivector annual_catch(1,Nyrs) // vector of total catches by year

 !!cout<< "Nages " << Nages <<endl;  
 !!cout<< "Nsizes " << Nsizes <<endl;
 !!cout<< "Nyrs " << Nyrs <<endl;  
 !!cout<< "nijk " << nijk <<endl;  
 !!cout<< "yjk " << yjk <<endl;
 
 //indices for non zero entries and ragged arrays
 init_int vec_length
 init_ivector k_vec(1,vec_length)
 init_ivector i_vec(1,vec_length)
 init_ivector j_vec(1,vec_length)
 init_ivector mini_vec(1,vec_length)
 init_ivector maxi_vec(1,vec_length)
 init_int vec_length2
 init_ivector k_vec2(1,vec_length2)
 init_ivector j_vec2(1,vec_length2)
 init_ivector mini_vec2(1,vec_length2)
 init_ivector maxi_vec2(1,vec_length2)


 !!cout<< "vec_length " << vec_length <<endl;
 !!cout<< "vec_length2 " << vec_length2 <<endl;
 
 int minPji
 !! minPji=1;
 !!cout<< "minPji " << minPji <<endl;

 int maxPji
 !! maxPji=Nages;
 !!cout<< "maxPji " << maxPji <<endl;


 init_ivector minind(minPji,maxPji) //min column indice for each row of Pji
 init_ivector maxind(minPji,maxPji) //max column indice for each row of Pji
 !!cout<< "minind " << minind <<endl;
 !!cout<< "maxind " << maxind <<endl;
 ivector rowlengths(1,Nages)
 !! rowlengths=maxind-minind+1;
 !!cout<< "rowlengths " << rowlengths <<endl;

PARAMETER_SECTION

 //Probability of age i  in year k
 init_matrix Pik_init(1,Nyrs,1,Nages-1)

 //Probability of size j given age i
 init_matrix Pji_init(minPji,maxPji,minind,maxind)

 matrix Pjiragged(minPji,maxPji,minind,maxind+1)
 
 !!cout<< "Pik_init" << Pik_init <<endl;
 !!cout<< "Pji_init" << Pji_init <<endl;

 sdreport_matrix Pik(1,Nyrs,1,Nages)
 sdreport_matrix Pji(1,Nages,1,Nsizes)
 sdreport_matrix CAA(1,Nyrs,1,Nages)
 
 objective_function_value nll // negative log likelihood


PROCEDURE_SECTION
 int indk;
 int indi;
 int indj;
 int indk2;
 int indj2;
 int indi1;
 int indi2;
 int indi11;
 int indi22;
 
  //filling Pji
 for(int i=1;i<=Nages;i++)
 {
  if(rowlengths(i)==1){
     Pjiragged.rowfill(i,Pji2qnum(Pji_init(i),minind[i]));
  } else {
     Pjiragged.rowfill(i,Pji2qtop(Pji_init(i),minind[i],maxind[i]));
  }
 }

 //filling Pik
 for(int k=1;k<=Nyrs;k++)
 {
   Pik.rowfill(k,Pik2p(Pik_init(k)));
 }


 //initialize variables
 nll.initialize();


 //likelihood equation 3 in Hoenig et al. 2002
 for(int l=1;l<=vec_length;l++)
 {
   indk=k_vec(l);
   indi=i_vec(l);
   indj=j_vec(l);
   indi11=mini_vec(l);
   indi22=maxi_vec(l);
   nll += log(Pik(indk,indi)*Pjiragged(indi,indj))*nijk(indk)(indi,indj);
 }

 for(int m=1;m<=vec_length2;m++)
 {
   indk2=k_vec2(m);
   indj2=j_vec2(m);
   indi1=mini_vec2(m);
   indi2=maxi_vec2(m);
   dvar_vector Pjicol(indi1,indi2);
   Pjicol.initialize();
   for (int n=indi1;n<=indi2;n++){
     Pjicol(n) = Pjiragged(n,indj2);
   }
   nll += log(Pjicol*Pik(indk2)(indi1,indi2))*yjk(indk2,indj2);
 }

 for(int k=1;k<=Nyrs;k++)
 {
   for(int i=1;i<=Nages;i++)
   {
     if(xik(k,i)==0){
     nll += 0;
     } else {
     nll += log(Pik(k,i))*xik(k,i);
     }
   }
 }
 
 nll *= -1;


//calculate CAA
  for(int k=1;k<=Nyrs;k++)
  {
   for(int i=1;i<=Nages;i++)
   {
    CAA(k,i)=Pik(k,i)*annual_catch(k);
   }
  }

//fill Pji full matrix
  for(int l=1;l<=vec_length;l++)
   {
    indi=i_vec(l);
    indj=j_vec(l);
    Pji(indi,indj)= Pjiragged(indi,indj);
   }

FUNCTION dvar_vector Pik2p(const dvar_vector& a) //use logit transform to get probability vector
  dvar_vector p(1,Nages);
  dvar_vector expa=mfexp(a);
  dvar_vector totalp=expa/(1+sum(expa)); 
  p(1,Nages-1)= totalp;
  p(Nages)=1-sum(p(1,Nages-1)); //last entry calculated as a fn of rest of the vector
  return p;

FUNCTION dvar_vector Pji2qtop(const dvar_vector& b, const int& bb,const int& bbb)
  dvar_vector q(bb,bbb+1);
  dvar_vector expb=mfexp(b);
  q(bb,bbb)=expb/(1+sum(expb));
  q(bbb+1)=1-sum(q(bb,bbb));
  return q;

FUNCTION dvar_vector Pji2qnum(const dvar_vector& c, const int& cc)
  int int0=cc;
  dvar_vector r(cc,cc+1);
  dvariable expc=mfexp(c[1]);
  r(cc)=expc/(1+expc);
  r(cc+1)=1-r(cc);
  return r;


RUNTIME_SECTION

 maximum_function_evaluations 2e5
 convergence_criteria 1e-6

REPORT_SECTION

 //save_gradients(gradients);

 report.precision(10);

 report << "# Nages" << endl;
 report << Nages << endl;
 report << " " << endl;

 report << "# Nsizes" << endl;
 report << Nsizes << endl;
 report << " " << endl;

 report << "# Nyrs" << endl;
 report << Nyrs << endl;
 report << " " << endl;
 
 report << "# nijk" << endl;
 report << nijk << endl;
 report << " " << endl;
 
 report << "# yjk" << endl;
 report << yjk << endl;
 report << " " << endl;
   
 report << "# xik" << endl;
 report << xik << endl;
 report << " " << endl;
 
 report << "# Pik" << endl;
 report << Pik << endl;
 report << " " << endl;
 
 report << "# Pji" << endl;
 report << Pji << endl;
 report << " " << endl;

 report << "# CAA" << endl;
 report << CAA << endl;
 report << " " << endl;
 
 report << "# nll" << endl;
 report << nll << endl;
 report << " " << endl;
  
 report << "# maximum gradient component" << endl;
 report << objective_function_value::pobjfun->gmax << endl;
 report << " " << endl;
