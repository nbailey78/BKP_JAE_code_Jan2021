new;

/*Bailey, N., Kapetanios, G. and Pesaran, H. P. (2021). Measurement of Factor Strength: Theory and Practice. 
                                                   Journal of Applied Econometrics, forthcoming.*/
/* Experiment 3B: Two unobserved factors - Gaussian errors (Bias, RMSE - two factors seq CSA) (d=0.25) */
/* Table 4 - main paper */

_dxmiss=0/0;
_optnorm=0; /* =0 normal errors; =1 t-distributed errors; =2 chi-squared errors */

repl=2000;

rndseed 12345;

grid_n=100|200|500|1000;
grid_t=60|120|200|500|1000;

p_val=0.1; /* p-value for the critical value c_p(n) */
a_size=0.1; /* p-value for JAE estimator */
p_test = 0.05; /* significance level of the test of H0: alpha=alpha0 */
delta =0.25; /* critical value exponent */

factor=2; /* =1, one factor; =2, two factors */
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


grid_alpha1=seqa(0.90,0.1,2).*.ones(3,1); /* strength of first factor */
grid_pr2={0.51,0.75,0.90};
grid_alpha2=ones(2,1).*.grid_pr2; /* stength of second factor */
grid_alpha=grid_alpha1~grid_alpha2;

bin=50;

all_o1={};
allbias_o1={};
allrmse_o1={};
all_o2={};
allbias_o2={};
allrmse_o2={};

nin=1;
do while nin<=rows(grid_n);
n=grid_n[nin];
tit=1;
do while tit<=rows(grid_t);
t=grid_t[tit];

sto_a1=zeros(repl,rows(grid_alpha));
sto_bias_o1=zeros(repl,rows(grid_alpha));
sto_mse_o1=zeros(repl,rows(grid_alpha));
sto_a2=zeros(repl,rows(grid_alpha));    
sto_bias_o2=zeros(repl,rows(grid_alpha));
sto_mse_o2=zeros(repl,rows(grid_alpha));    


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
  zita=rndMVn(t+bin, m_f, cov_f);
  /*zita=sig_f*rndn(t+bin,factor);*/  
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
a02=grid_alpha[m,2];
a0=a01~a02;

/* Construction of loadings */
gamm=zeros(factor,n);
gamm[1,1:floor(n^a01)]=gamm1[1,1:floor(n^a01)];
gamm[2,1:floor(n^a02)]=gamm2[1,1:floor(n^a02)]; 
gamm_adj=gamm[2,.];
iy=rndu(cols(gamm_adj),1);
h_mix=sorthc(gamm_adj'~iy,2);
gamm[2,.]=(h_mix[.,1])';
  
/* Dataset */
x=ones(t,n)*c_v+f*gamm+eps;
 
 /* Estimation of alpha */
 /* Sequential CSA */
{a1,a2}=seq(x,p_val);

/* Save */
  sto_a1[o,m]=a1;
  sto_bias_o1[o,m]=a1-a01;
  sto_mse_o1[o,m]=(a1-a01)^2;  
  sto_a2[o,m]=a2;
  sto_bias_o2[o,m]=a2-a02;
  sto_mse_o2[o,m]=(a2-a02)^2;  
endfor;
o=o+1;
endo;
m_o1=(sumc(sto_a1)'/repl);
m_bias_o1=(sumc(sto_bias_o1)'/repl)*100;
m_mse_o1=(sqrt(sumc(sto_mse_o1)'/repl))*100;
m_o2=(sumc(sto_a2)'/repl);
m_bias_o2=(sumc(sto_bias_o2)'/repl)*100;
m_mse_o2=(sqrt(sumc(sto_mse_o2)'/repl))*100;

all_o1=all_o1|(m_o1);
allbias_o1=allbias_o1|(m_bias_o1);
allrmse_o1=allrmse_o1|(m_mse_o1);
all_o2=all_o2|(m_o2);
allbias_o2=allbias_o2|(m_bias_o2);
allrmse_o2=allrmse_o2|(m_mse_o2);

tit=tit+1;
endo;
nin=nin+1;
endo;

"all_o1";all_o1;
"allbias_o1";allbias_o1;
"allrmse_o1";allrmse_o1;
"all_o2";all_o2;
"allbias_o2";allbias_o2;
"allrmse_o2";allrmse_o2;


/* *********************PROCEDURES********************* */

proc (2)=seq(x,p_val);
local n,csa,pi_hat,a1,ucsa,csau,pi_res,a2;
n=cols(x);
csa=sumc((rndu(n,1)+1).*(x'));
pi_hat=meanc(testm(x,csa,p_val,delta));    
if pi_hat==0;
a1=0;  
else;
a1=1+(ln(pi_hat)/ln(cols(x)));
endif;
ucsa=x-csa*(x/csa);
csau=sumc((rndu(n,1)+2).*(ucsa'));
pi_res=meanc(testm(ucsa,csau,p_val,delta));
if pi_res==0;
a2=0;
else;    
a2=1+(ln(meanc(testm(ucsa,csau,p_val,delta)))/ln(cols(x)));
endif;
retp(a1,a2);
endp;

proc testm(y,x,p,delta);
local xx,b,pv,alpha,t,n;
y=y-meanc(y)';
x=x-meanc(x)';
t=rows(x);
n=cols(y);
xx=sumc(x.*x);
b=sumc(y.*x)./xx;
pv=cdfn(b./sqrt(meanc((y-x.*(b'))^2)./xx));
alpha=p/(n^delta);
retp((pv.<alpha).or((1-pv).<alpha));
endp;
