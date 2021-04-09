new;

/*Bailey, N., Kapetanios, G. and Pesaran, H. P. (2021). Measurement of Factor Strength: Theory and Practice. 
                                                   Journal of Applied Econometrics, forthcoming.*/
/* Experiment 1A: Observed single factor - Gaussian errors (Bias, RMSE, Size) */
/* Table S1a - appendix */

_dxmiss=0/0;
_optnorm=0; /* =0 normal errors; =1 t-distributed errors; =2 chi-squared errors */

repl=2000;

rndseed 12345;

grid_n=100|200|500|1000;
grid_t=60|120|200|500|1000;

p_val=0.1; /* p-value for the critical value c_p(n) */
p_test = 0.05; /* significance level of the test of H0: alpha=alpha0 */
delta =0.25; /* critical value exponent */

factor=1; /* =1, one factor; =2, two factors */
rho_f=0.5; /* AR strength of factors */
sig_f=1; /* variance of factor */
cov_f = { 1 0.0, 0.0 1 }; /* Covarianve structure of factors: 0, 0.25 */
m_f = { 0, 0 }; /* Mean of factors */
mu_c=0; /* mean of constant */
mu1=sqrt(0.75); /* mean of 1st factor's loadings */
mu2=sqrt(0.75);  /* mean of 1st factor's loadings */
mu1_b=0.2; /* Bound of first factor loadings */
mu2_b=0.2; /* Bound of second factor loadings */
sc_sige=1; /* or 3/4: scale of variance of errors */
v_chi=2; /* Chi-sqr degrees of freedom for variance of errors */
v_chi1=8; /* Chi-sqr degrees of freedom for t-distributed errors */

grid_alpha1=seqa(0.7,0.05,7).*.ones(7,1); /* strength of first factor */
grid_alpha2=ones(7,1).*.seqa(0.7,0.05,7); /* stength of second factor */
if factor==1;
grid_alpha=seqa(0.7,0.05,7);
elseif factor==2; 
grid_alpha=grid_alpha1~grid_alpha2;
endif;

bin=50;

allbias_o1={};
allrmse_o1={};
allsize_fo1={};
allpowerp_fo1={};
allpowerm_fo1={};

allbias_o2={};
allrmse_o2={};
allsize_fo2={};
allpowerp_fo2={};
allpowerm_fo2={};

nin=1;
do while nin<=rows(grid_n);
n=grid_n[nin];
tit=1;
do while tit<=rows(grid_t);
t=grid_t[tit];

sto_bias_o1=zeros(repl,rows(grid_alpha));
sto_mse_o1=zeros(repl,rows(grid_alpha));    
sto_size_fo1=zeros(repl,rows(grid_alpha));
sto_powerp_fo1=zeros(repl,rows(grid_alpha));
sto_powerm_fo1=zeros(repl,rows(grid_alpha));    

sto_bias_o2=zeros(repl,rows(grid_alpha));
sto_mse_o2=zeros(repl,rows(grid_alpha));   
sto_size_fo2=zeros(repl,rows(grid_alpha));
sto_powerp_fo2=zeros(repl,rows(grid_alpha));
sto_powerm_fo2=zeros(repl,rows(grid_alpha));    
   
o=1;
do while  o<=repl; /*loop of replications*/
  
/* Construction of error term */
if _optnorm==0; /* Normal errors */ 
e=rndn(t+bin,n);
chisq=1/3*(1+(sumc(rndn(n,v_chi)'.^2)));
chisq=diagrv(eye(n),chisq);    
eps=e*(chisq.^0.5)';    
elseif _optnorm==1; /* t-distributed errors */
e=rndn(t+bin,n);    
chisq=1/3*(1+(sumc(rndn(n,v_chi)'.^2)));
chisq=diagrv(eye(n),chisq); 
e1=e*(chisq.^0.5)';      
chisq_t=(sumc(rndn(t+bin,v_chi1)'.^2))';
chisq1=ones(n,t+bin).*chisq_t;
eps=(((v_chi1-2)./chisq1).^0.5)'.*e1;
elseif _optnorm==2; /* chi-squared distributed errors */
e=(rndn(t+bin,n).^2+rndn(t+bin,n).^2); /* 2 degrees of freedom */
chisq=1/3*(1+(sumc(rndn(n,v_chi)'.^2)));
chisq=diagrv(eye(n),chisq);    
eps=(e-2)*(chisq.^0.5)/2';      
endif;    
eps=eps[bin+1:t+bin,.];

/* Construction of factors */
if factor==1;
  zita=sig_f*rndn(t+bin,factor);
elseif factor==2;
  /*zita=rndMVn(t+bin, m_f, cov_f);*/
  zita=sig_f*rndn(t+bin,factor);  
endif;
f=zeros(t+bin,factor);
if rho_f==1;
 for p(2,t+bin,1);
  f[p,.]=rho_f*f[p-1,.]+zita[p,.];
 endfor;    
else;
 for p(2,t+bin,1);
  f[p,.]=rho_f*f[p-1,.]+sqrt((1-rho_f^2))*zita[p,.];
 endfor;
endif;
f=f[bin+1:t+bin,.];

/* Constant */
const=mu_c+rndn(n,1);
c_v=diagrv(eye(n),const);

/* Construction of V_i's */
gamm1=(mu1-mu1_b)+(2*mu1_b)*rndu(1,n); 
if factor==2;
gamm2=(mu2-mu2_b)+(2*mu2_b)*rndu(1,n);
endif;

for m(1,rows(grid_alpha),1);  /* loop of strength of factors */
  n~t~o~grid_alpha[m,.];

/* Strength of factors */
a01=grid_alpha[m,1];
a1p1=a01+0.05;   
if a1p1>0.96;
 a1p1=0.98;  
endif;
a1m1=a01-0.05;

if factor==1;
 a0=a01;
 a1p=a1p1;
 a1m=a1m1;   
elseif factor==2;
 a02=grid_alpha[m,2];
 a1p2=a02+0.05;      
  if a1p2>0.96;
   a1p2=0.98;  
  endif;
 a1m2=a02-0.05;
 a0=a01~a02;
 a1p=a1p1~a1p2;
 a1m=a1m1~a1m2;      
endif;

/* Construction of loadings */
gamm=zeros(factor,n);
gamm1p=zeros(factor,n);
gamm1m=zeros(factor,n);
 if factor==1;
  gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];
  gamm1p[1,1:floor(n^a1p1)]=gamm1[1,1:floor(n^a1p1)];
  gamm1m[1,1:floor(n^a1m1)]=gamm1[1,1:floor(n^a1m1)];   
 elseif factor==2;
  gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];
  gamm1p[1,1:floor(n^a1p1)]=gamm1[1,1:floor(n^a1p1)];
  gamm1m[1,1:floor(n^a1m1)]=gamm1[1,1:floor(n^a1m1)];      
  gamm[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)]; 
  gamm1p[2,1:floor(n^a1p2)]=gamm2[1,1:floor(n^a1p2)];
  gamm1m[2,1:floor(n^a1m2)]=gamm2[1,1:floor(n^a1m2)]; 
  gamm_adj=gamm[2,.]|gamm1p[2,.]|gamm1m[2,.];
  iy=rndu(cols(gamm_adj),1);
  h_mix=sorthc(gamm_adj'~iy,4);
  gamm[2,.]=(h_mix[.,1])';
  gamm1p[2,.]=(h_mix[.,2])';
  gamm1m[2,.]=(h_mix[.,3])';
 endif;
 
/* Dataset */
x=ones(t,n)*c_v+f*gamm+eps;
x1p=ones(t,n)*c_v+f*gamm1p+eps;
x1m=ones(t,n)*c_v+f*gamm1m+eps; 
 
 /* Estimation of alpha */
 /* Observed factor */
{a_hato,size_fo,rej,rej}=a_d(x,f,p_val,delta,p_test,a0,a1p,a1m,0);
/*{rej,rej,powerp_fo,rej}=a_d(x1p,f,p_val,delta,p_test,a0,a1p,a1m,0); 
{rej,rej,rej,powerm_fo}=a_d(x1m,f,p_val,delta,p_test,a0,a1p,a1m,0);*/ 

/* Save */
if factor==1;
  sto_bias_o1[o,m]=a_hato-a0;
  sto_mse_o1[o,m]=(a_hato-a0)^2;
  sto_size_fo1[o,m]=size_fo;
  /*sto_powerp_fo1[o,m]=powerp_fo;
  sto_powerm_fo1[o,m]=powerm_fo; */   
 elseif factor==2;
  sto_bias_o1[o,m]=a_hato[1,1]-a01;
  sto_mse_o1[o,m]=(a_hato[1,1]-a01)^2;
  sto_bias_o2[o,m]=a_hato[1,2]-a02;
  sto_mse_o2[o,m]=(a_hato[1,2]-a02)^2;
  sto_size_fo1[o,m]=size_fo[1,1];
  /*sto_powerp_fo1[o,m]=powerp_fo[1,1];
  sto_powerm_fo1[o,m]=powerm_fo[1,1];*/
  sto_size_fo2[o,m]=size_fo[1,2];
  /*sto_powerp_fo2[o,m]=powerp_fo[1,2];
  sto_powerm_fo2[o,]=powerm_fo[1,2];*/
endif;
endfor;
o=o+1;
endo;

m_bias_o1=(sumc(sto_bias_o1)'/repl);
m_mse_o1=(sqrt(sumc(sto_mse_o1)'/repl));
m_size_fo1=(sumc(sto_size_fo1)'/repl);
/*m_powerp_fo1=(sumc(sto_powerp_fo1)'/repl);
m_powerm_fo1=(sumc(sto_powerm_fo1)'/repl);*/

m_bias_o2=(sumc(sto_bias_o2)'/repl);
m_mse_o2=(sqrt(sumc(sto_mse_o2)'/repl));
m_size_fo2=(sumc(sto_size_fo2)'/repl);
/*m_powerp_fo2=(sumc(sto_powerp_fo2)'/repl);
m_powerm_fo2=(sumc(sto_powerm_fo2)'/repl);*/

allbias_o1=allbias_o1|(m_bias_o1);
allrmse_o1=allrmse_o1|(m_mse_o1);
allsize_fo1=allsize_fo1|(m_size_fo1);
/*allpowerp_fo1=allpowerp_fo1|(m_powerp_fo1);
allpowerm_fo1=allpowerm_fo1|(m_powerm_fo1);*/

allbias_o2=allbias_o2|(m_bias_o2);
allrmse_o2=allrmse_o2|(m_mse_o2);
allsize_fo2=allsize_fo2|(m_size_fo2);
/*allpowerp_fo2=allpowerp_fo2|(m_powerp_fo2);
allpowerm_fo2=allpowerm_fo2|(m_powerm_fo2);*/

tit=tit+1;
endo;
nin=nin+1;
endo;


"allbias_o1";allbias_o1;
"allrmse_o1";allrmse_o1;
"allsize_fo1";allsize_fo1;
/*"allpowerp_fo1";allpowerp_fo1;
"allpowerm_fo1";allpowerm_fo1;*/
"";"";
"allbias_o2";allbias_o2;
"allrmse_o2";allrmse_o2;
"allsize_fo2";allsize_fo2;
/*"allpowerp_fo2";allpowerp_fo2;
"allpowerm_fo2";allpowerm_fo2;*/


/* *********************PROCEDURES********************* */

proc (4)=a_d(x,f,p_val,delta,p_test,a0,a1p,a1m,opt); /* Dataset (TxN); true factors (Txr); size of MT (eg 0.10); MT adjustment (eg 1/3); size of a_hat testing; alpha under: null, alt+, alt-, observed(=0)/ unobserved (=1) factor */
local n,t,x_ave,rhsv,coef,res,sig_vec,t_stat,cv,di,d_bar,a_hat,ze_nf,psi_nf,z_af,z_a1pf,z_a1mf,ze_ninf,psi_ninf,ze_ninf1p,psi_ninf1p,ze_ninf1m,psi_ninf1m,z_ainf,z_a1pinf,z_a1minf,size_f,power_pf,power_mf,
    size_inf,power_pinf,power_minf,crit,crit1,a_adj,r,level,fi,ur1,ur2,ur,csi,csi2,thetavec,iu,ir,thetavecdm,thetau,thetau2,thetamean;
n=cols(x);
t=rows(x);    
r=400;
crit=cdfni(1-(p_test/2));
level=p_test/(8*ceil((n/100)^(1/4)));      
/* OLS regression */
if opt==0; /* use observed f */  
 x_ave=f;
 elseif opt==1; /* use cross-sectional averages */
 x_ave=meanc(x');    
endif;
rhsv=ones(t,1)~x_ave;
coef=inv(rhsv'*rhsv)*rhsv'*x; /* OLS estimate from regression of x on a constant and the factor(s) */
res=x-rhsv*coef;
sig_vec=res'*res/(t-cols(rhsv));
t_stat=coef[2:rows(coef),.]./(sqrt(diag(sig_vec))'.*sqrt(diag(inv(x_ave'*x_ave))));
cv=cdfni(1-(p_val/(2*(n^delta))));
di=(abs(t_stat).>cv)';
d_bar=sumc(di)';
a_hat=zeros(1,cols(d_bar));
z_af=zeros(1,cols(d_bar));
z_a1pf=zeros(1,cols(d_bar));
z_a1mf=zeros(1,cols(d_bar));
size_f=zeros(1,cols(d_bar));
power_pf=zeros(1,cols(d_bar));
power_mf=zeros(1,cols(d_bar));  
fi=zeros(1,cols(d_bar));
ur1=-exp(1/0.03);
ur2=exp(1/0.03);
for j(1,cols(d_bar),1);
  if d_bar[1,j]==0;
   a_hat[1,j]=0;
  else;
   a_hat[1,j]=ln(d_bar[1,j])/ln(n);
  endif;
 /* feasible */
 ze_nf=(p_val*(n-n^(a_hat[1,j])))/(ln(n)*n^(delta+a_hat[1,j]));
 psi_nf=p_val*(n-n^(a_hat[1,j]))*n^(-delta-2*a_hat[1,j])*(1-(p_val/n^(delta)));
 if a_hat[1,j]==1;
  z_af[1,j]=999;
  z_a1pf[1,j]=999;
  z_a1mf[1,j]=999;   
 else;  
  z_af[1,j]=(ln(n)*(a_hat[1,j]-a0[1,j]-ze_nf))/sqrt(psi_nf);
  z_a1pf[1,j]=(ln(n)*(a_hat[1,j]-a0[1,j]-ze_nf))/sqrt(psi_nf);
  z_a1mf[1,j]=(ln(n)*(a_hat[1,j]-a0[1,j]-ze_nf))/sqrt(psi_nf);   
 endif; 
 size_f[1,j]=(abs(z_af[1,j])>crit);
 power_pf[1,j]=(abs(z_a1pf[1,j])>crit);
 power_mf[1,j]=(abs(z_a1mf[1,j])>crit);  
 if a0[1,j]>0.999;
  if a0[1,j]==a_hat[1,j];
   fi[1,j]=exp(t);
  else;
   fi[1,j]=(1/(abs(a0[1,j]-a_hat[1,j])));   
  endif;   
  csi=rndn(r,1);
  csi2=(fi[1,j])*csi;
  ur=ur1|ur2;   
  thetavec=zeros(r,rows(ur));
  iu=1;
   do while iu<=rows(ur);
   ir=1;
    do while ir<=r;
     if csi2[ir,.]<ur[iu,.];
      thetavec[ir,iu]=1;
     endif;                
    ir=ir+1;
   endo;
  iu=iu+1;
  endo;
  thetavecdm=thetavec-(0.5)*ones(r,rows(ur));
  thetau=(2/sqrt(r))*sumc(thetavecdm); 
  thetau2=thetau^2;
  thetamean=sumc(thetau2/rows(ur)); 
  size_f[1,j]=(thetamean>cdfchii(1-level,1));  
 endif;
endfor;
retp(a_hat,size_f,power_pf,power_mf);  
endp;
