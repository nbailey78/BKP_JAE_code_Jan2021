new;
cls;

p_val=0.1;
p_test=0.05;
delta=1/4;
p_v=12;

a_size=0.05;
a_size1=0.10;
delta_thr=1/2;
delta2_thr=1/4;

P=340; /* no. of monthly S&P compositions - 340 in total (until Dec 2017) */

all_amkt14={};
all_amktcip14={};
all_amktcim14={};
all_ahat14={};
all_adj14={};
all_acip14={};
all_acim14={};

sto_n_or=zeros(P,1);
sto_t_or=zeros(P,1);
sto_n=zeros(P,1);
sto_t=zeros(P,1);

for i(1,P,1);

/* PART I: UPLOAD AND CLEAN THE DATA */
/* Load Data */
b="SP";
e=".xls";
f_number=ftos(i,"%lf",1,0);
data=b$+f_number$+e;

/* Choice of starting point: 84 for 120M hist, 144 for 60M hist */
x1 = xlsreadm(data,"b84:fy542",1,0);
x2 = xlsreadm(data,"b84:fy542",2,0);
x3 = xlsreadm(data,"b84:ek542",3,0);
rf=xlsreadm("RF.xls","b84:b542",1,0);
mkt_rf=xlsreadm("Mkt_RF2.xls","b84:eq542",1,0);/* 147 factors: end in Dec 2017 (point 542) */   
mkt_rf=mkt_rf;    
x_loop=x1~x2~x3;

t_adj=120; /* Choice of history in no. of months: 120, 60 */

x_or=x_loop[i:i+t_adj-1,.];
rf_loop=rf[i:i+t_adj-1,1];
mkt_ff=mkt_rf[i:i+t_adj-1,.];

t_or=rows(x_or);
n_or=cols(x_or);
sto_n_or[i,1]=n_or;
sto_t_or[i,1]=t_or;

/* Delete stocks without required history */
x_trans=x_or';
e0=(x_trans[.,1] .eq 0);
x_adj=delif(x_trans,e0);
/* Delete stocks without current data */
e1=(x_adj[.,t_adj] .eq 0);
x_adj1=delif(x_adj,e1);
/* Delete stocks without 5 successive periods of data */
for o(1,rows(x_adj1),1);
for j(1,cols(x_adj1)-4,1);
    if x_adj1[o,j]==0 and x_adj1[o,j+1]==0 and x_adj1[o,j+2]==0 and x_adj1[o,j+3]==0 and x_adj1[o,j+4]==0;
     x_adj1[o,1]=-999;
    endif;
endfor;
endfor;
e2=(x_adj1[.,1] .eq -999);
x_adj2=delif(x_adj1,e2);

/* Final dataset: x_prel=x_r if simple returns, x=x_r-rf_loop if excess returns, x_stand=x_stand_prel-mkt_stand1*coef if residual returns */
x_r=x_adj2';
/* Excess returns */
x=x_r-rf_loop;
m_x=meanc(x)';
std_x=stdc(x)';

t=rows(x);
n=cols(x);
sto_n[i,1]=n;
sto_t[i,1]=t;

/* Residuals - excluding mkt */
interc=ones(t,1);
coef=inv(interc'*interc)*interc'*x;
x_prel=x-interc*coef; 

/* PART II: ESTIMATION OF FACTOR STRENGTH */ 
/* Estimation of rolling factor strengths for 146 factors */
a_mkt14=zeros(1,cols(mkt_ff)-1);
a_mktcip14=zeros(1,cols(mkt_ff)-1);
a_mktcim14=zeros(1,cols(mkt_ff)-1);
a_hat14=zeros(1,cols(mkt_ff)-1);
a_cip14=zeros(1,cols(mkt_ff)-1);
a_cim14=zeros(1,cols(mkt_ff)-1);
for j(2,cols(mkt_ff),1);
  ff=mkt_ff[.,1]~mkt_ff[.,j];
  {a_hat14p,rej,ci_p14,ci_m14,mR2_14}=a_demp(x,ff,p_val,delta,p_test);
  a_mkt14[1,j-1]=a_hat14p[1,1];
  a_mktcip14[1,j-1]=ci_p14[1,1];
  a_mktcim14[1,j-1]=ci_m14[1,1];  
  a_hat14[1,j-1]=a_hat14p[1,2];
  a_cip14[1,j-1]=ci_p14[1,2];
  a_cim14[1,j-1]=ci_m14[1,2];    
endfor;

/* Regression with just the mkt factor or unobserved factor */
  /*ff=mkt_ff[.,1];*/ /* Mkt factor */
  /*ff=meanc(x');*/ /* Unobserved factos (CSA) */
  /*{a_hat14p,rej,ci_p14,ci_m14,mR2_14}=a_demp(x,ff,p_val,delta,p_test);
  a_mkt14=a_hat14p[1,1];
  a_mktcip14=ci_p14[1,1];
  a_mktcim14=ci_m14[1,1]; */


/* Save */
all_amkt14=all_amkt14|a_mkt14;
all_amktcip14=all_amktcip14|a_mktcip14;
all_amktcim14=all_amktcim14|a_mktcim14;
all_ahat14=all_ahat14|a_hat14;
all_acip14=all_acip14|a_cip14;
all_acim14=all_acim14|a_cim14;

endfor;

"n";sto_n_or~sto_n;
"t";sto_t_or~sto_t;
"sto_amkt14";all_amkt14;
"sto_amktcim14";all_amktcim14;
"sto_amktcip14";all_amktcip14;
"sto_ahat14";all_ahat14;
"sto_acim14";all_acim14;
"sto_acip14";all_acip14;

/* *********************PROCEDURES********************* */

proc (5)=a_demp(x,x_ave,p_val,delta,p_test); /* Dataset (TxN); factor; size of MT (eg 0.10); MT adjustment (eg 1/3); size of a_hat testing */
local t,n,rhsv,coef,res,m_x,x_dm,R2,mR2,sig_vec,t_stat,cv,di,d_bar,a_hat,a_adj,ze_nf,psi_nf,ci_a,ci_p,ci_m;
t=rows(x);
n=cols(x);
rhsv=ones(t,1)~x_ave;
coef=inv(rhsv'*rhsv)*rhsv'*x; /* OLS estimate from regression of x on a constant and the factor */
res=x-rhsv*coef;
m_x=ones(t,n)*diagrv(eye(n),meanc(x)); 
x_dm=x-m_x;
R2=1-diag(res'*res)./diag(x_dm'*x_dm);
mR2=meanc(R2);    
sig_vec=res'*res/(t-cols(rhsv));
t_stat=coef[2:rows(coef),.]./(sqrt(diag(sig_vec))'.*sqrt(diag(inv(x_ave'*x_ave))));
cv=cdfni(1-(p_val/(2*(n^delta))));
di=(abs(t_stat).>cv)';
d_bar=sumc(di)';
a_hat=zeros(1,cols(d_bar));
a_adj=zeros(1,cols(d_bar));
ze_nf=zeros(1,cols(d_bar));
psi_nf=zeros(1,cols(d_bar));
ci_a=zeros(1,cols(d_bar));
ci_p=zeros(1,cols(d_bar));
ci_m=zeros(1,cols(d_bar));
for j(1,cols(d_bar),1);
if d_bar[1,j]==0;
 a_hat[1,j]=0;
else;
 a_hat[1,j]=ln(d_bar[1,j])/ln(n);
endif;
ze_nf[1,j]=(p_val*(n-n^(a_hat[1,j])))/(ln(n)*n^(delta+a_hat[1,j]));
psi_nf[1,j]=p_val*(n-n^(a_hat[1,j]))*n^(-delta-2*a_hat[1,j])*(1-(p_val/n^(delta)));
a_adj[1,j]=a_hat[1,j]-ze_nf[1,j];
/* 90% CI: 1.65, 95% CI: 1.96 */
ci_a[1,j]=1.65*(sqrt(psi_nf[1,j])/ln(n));
ci_p[1,j]=a_hat[1,j]+ci_a[1,j];
ci_m[1,j]=a_hat[1,j]-ci_a[1,j];
endfor;
retp(a_hat,a_adj,ci_p,ci_m,mR2);  
endp;

