new;

/*Bailey, N., Kapetanios, G. and Pesaran, H. P. (2021). Measurement of Factor Strength: Theory and Practice. 
                                                   Journal of Applied Econometrics, forthcoming.*/
/* Experiment 1A: Observed single factor - Gaussian errors (Power functions) */
/* Figure S1a - appendix */

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
cov_f = { 1 0.0, 0.0 1 }; /* Covarianve structure of factors */
m_f = { 0, 0 }; /* Mean of factors */
mu_c=0; /* mean of constant */
mu1=sqrt(0.75); /* mean of 1st factor's loadings */
mu2=sqrt(0.75);  /* mean of 1st factor's loadings */
mu1_b=0.2; /* Bound of first factor loadings */
mu2_b=0.2; /* Bound of second factor loadings */
sc_sige=1; /* or 3/4: scale of variance of errors */
v_chi=2; /* Chi-sqr degrees of freedom for variance of errors */
v_chi1=8; /* Chi-sqr degrees of freedom for t-distributed errors */
bin=50;

grid_alpha1=seqa(0.7,0.05,7).*.ones(7,1); /* strength of first factor */
grid_alpha2=ones(7,1).*.seqa(0.7,0.05,7); /* stength of second factor */
if factor==1;
grid_alpha=seqa(0.7,0.05,7);
elseif factor==2; 
grid_alpha=grid_alpha1~grid_alpha2;
endif;
ap1=seqa(-0.05,0.005,21); /* grid for power - f1 */
ap2=seqa(-0.05,0.005,21); /* grid for power - f2 */   

allpowerp_fo1={};
allpowerp_fo2={};

powerp1=zeros(1,rows(ap1));
powerp2=zeros(1,rows(ap2));

nin=1;
do while nin<=rows(grid_n);
n=grid_n[nin];
tit=1;
do while tit<=rows(grid_t);
t=grid_t[tit];

m_power_o1={};
m_power_o2={};

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

m_power_ra1={};
m_power_ra2={};
power_ra1={};
power_ra2={};

 for m(1,rows(grid_alpha),1);  /* loop of strength of factors */
  n~t~o~grid_alpha[m,.];

 /* Strength of factors */
 a01=grid_alpha[m,1];
 if factor==1;
  a0=a01;
  ap=a0+ap1;
 elseif factor==2;
  a02=grid_alpha[m,2];
  a0=a01~a02;
  ap=(a01+ap1)~(a02+ap2);
 endif;

 /* Construction of loadings */
 gamm=zeros(factor,n);
 gammp=zeros(factor,n);
 for r(1,rows(ap),1);
  if factor==1;
   if a01>0.999999999999999 and ap[r,1]<1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];   
    gammp[1,1:floor(n^ap[r,1])]=gamm1[1,1:floor(n^ap[r,1])];
   elseif a01>0.999999999999999 and ap[r,1]>=1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];   
    gammp[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];     
   elseif a01<1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];
    gammp[1,1:floor(n^ap[r,1])]=gamm1[1,1:floor(n^ap[r,1])];
   endif;
 elseif factor==2;
  if a01>0.999999999999999 and ap[r,1]<1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];  
    gammp[1,1:floor(n^ap[r,1])]=gamm1[1,1:floor(n^ap[r,1])];
  elseif a01>0.999999999999999 and ap[r,1]>=1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];   
    gammp[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];     
   elseif a01<1;
    gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];
    gammp[1,1:floor(n^ap[r,1])]=gamm1[1,1:floor(n^ap[r,1])];
   endif;
  if a02>0.999999999999999 and ap[r,2]<1;
   gamm[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)]; 
   gammp[2,1:floor(n^ap[r,2])]=gamm2[1,1:floor(n^ap[r,2])];
  elseif a02>0.999999999999999 and ap[r,2]>=1;
   gamm[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)]; 
   gammp[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)];   
  elseif a02<1;
   gamm[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)]; 
   gammp[2,1:floor(n^ap[r,2])]=gamm2[1,1:floor(n^ap[r,2])];
  endif;
  gamm_adj=gamm[2,.]|gammp[2,.];
  iy=rndu(cols(gamm_adj),1);
  h_mix=sorthc(gamm_adj'~iy,3);
  gamm[2,.]=(h_mix[.,1])';
  gammp[2,.]=(h_mix[.,2])';
 endif;
 
 /* Dataset */
 x0=ones(t,n)*c_v+f*gamm+eps;
 xp=ones(t,n)*c_v+f*gammp+eps;
 
 /* Estimation of alpha */
 /* Observed factor */
 {rej,rej,power_fo}=a_d(x0,xp,f,p_val,delta,p_test,a0,ap[r,1],0); 

 if factor==1;
   powerp1[1,r]=power_fo[1,1]; 
  elseif factor==2;
   powerp1[1,r]=power_fo[1,1]; 
   powerp2[1,r]=power_fo[1,2];    
  endif;
 endfor;

/* Save */
  if factor==1;
    power_ra1=power_ra1|(powerp1);
   elseif factor==2;
    power_ra1=power_ra1|(powerp1);
    power_ra2=power_ra2|(powerp2);
  endif;
 endfor;
if factor==1;
 m_powerp1=reshape(power_ra1,1,rows(grid_alpha)*rows(ap));     
 m_power_o1=m_power_o1|(m_powerp1);
elseif factor==2; 
 m_powerp1=reshape(power_ra1,1,rows(grid_alpha)*rows(ap));  
 m_power_o1=m_power_o1|(m_powerp1);
 m_powerp2=reshape(power_ra2,1,rows(grid_alpha)*rows(ap));  
 m_power_o2=m_power_o2|(m_powerp2);
endif;
o=o+1;
endo;
if factor==1;
  m_powerp_fo1=sumc(m_power_o1)'/repl; 
  allpowerp_fo1=allpowerp_fo1|(m_powerp_fo1);  
elseif factor==2; 
  m_powerp_fo1=sumc(m_power_o1)'/repl;   
  m_powerp_fo2=sumc(m_power_o2)'/repl; 
  allpowerp_fo1=allpowerp_fo1|(m_powerp_fo1);
  allpowerp_fo2=allpowerp_fo2|(m_powerp_fo2);  
endif;   

tit=tit+1;
endo;
nin=nin+1;
endo;

"allpowerp_fo1";allpowerp_fo1;
"";"";
"allpowerp_fo2";allpowerp_fo2;


/* *********************PROCEDURES********************* */

proc (3)=a_d(y,x,f,p_val,delta,p_test,a0,ap,opt); /* Dataset (TxN); true factors (Txr); size of MT (eg 0.10); MT adjustment (eg 1/3); size of a_hat testing; alpha under: null, alt+, alt-, observed(=0)/ unobserved (=1) factor */
local n,t,x_ave,rhsv,coef,res,sig_vec,t_stat,cv,di,d_bar,a_hat,coef0,res0,sig_vec0,t_stat0,di0,d_bar0,a_hat0,ze_nf,psi_nf,z_af,z_a1pf,z_a1mf,size_f,power_f,power_mf,
    crit,crit1,a_adj,r,level,fi,ur1,ur2,ur,csi,csi2,thetavec,iu,ir,thetavecdm,thetau,thetau2,thetamean;
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
coef=inv(rhsv'*rhsv)*rhsv'*x; /* OLS estimate from regression of x under the alternative on a constant and the factor(s) */
res=x-rhsv*coef;
sig_vec=res'*res/(t-cols(rhsv));
t_stat=coef[2:rows(coef),.]./(sqrt(diag(sig_vec))'.*sqrt(diag(inv(x_ave'*x_ave))));
cv=cdfni(1-(p_val/(2*(n^delta))));
di=(abs(t_stat).>cv)';
d_bar=sumc(di)';

coef0=inv(rhsv'*rhsv)*rhsv'*y; /* OLS estimate from regression of x under the null on a constant and the factor(s) */
res0=y-rhsv*coef0;
sig_vec0=res0'*res0/(t-cols(rhsv));
t_stat0=coef0[2:rows(coef0),.]./(sqrt(diag(sig_vec0))'.*sqrt(diag(inv(x_ave'*x_ave))));
di0=(abs(t_stat0).>cv)';
d_bar0=sumc(di0)';

a_hat=zeros(1,cols(d_bar));
a_hat0=zeros(1,cols(d_bar0));
z_af=zeros(1,cols(d_bar));
z_a1pf=zeros(1,cols(d_bar));
size_f=zeros(1,cols(d_bar));
power_f=zeros(1,cols(d_bar));
fi=zeros(1,cols(d_bar));
ur1=-exp(1/0.03);
ur2=exp(1/0.03);
for j(1,cols(d_bar),1);
  if d_bar[1,j]==0;
   a_hat[1,j]=0;
  else;
   a_hat[1,j]=ln(d_bar[1,j])/ln(n);
  endif;
  if d_bar0[1,j]==0;
   a_hat0[1,j]=0;
  else;
   a_hat0[1,j]=ln(d_bar0[1,j])/ln(n);
  endif;
 /* feasible */
 ze_nf=(p_val*(n-n^(a_hat[1,j])))/(ln(n)*n^(delta+a_hat[1,j]));
 psi_nf=p_val*(n-n^(a_hat[1,j]))*n^(-delta-2*a_hat[1,j])*(1-(p_val/n^(delta)));
 if a_hat[1,j]==1;
  z_af[1,j]=999;
  z_a1pf[1,j]=999;
 else; 
  z_af[1,j]=(ln(n)*(a_hat[1,j]-a0[1,j]-ze_nf))/sqrt(psi_nf);
  z_a1pf[1,j]=(ln(n)*(a_hat[1,j]-a0[1,j]-ze_nf))/sqrt(psi_nf);
 endif;
 crit=cdfni(1-(p_test/2));
 size_f[1,j]=(abs(z_af[1,j])>crit);
 power_f[1,j]=(abs(z_a1pf[1,j])>crit);
 if a0[1,j]>0.999;
  if ap[1,j]==a_hat0[1,j];
   fi[1,j]=exp(t);
  else;
   fi[1,j]=(1/(abs(ap[1,j]-a_hat0[1,j])));   
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
  power_f[1,j]=(thetamean>cdfchii(1-level,1));  
 endif;
endfor;
retp(a_hat,size_f,power_f);  
endp;
