new;
cls;

/*Bailey, N., Kapetanios, G. and Pesaran, H. P. (2021). Measurement of Factor Strength: Theory and Practice. 
                                                   Journal of Applied Econometrics, forthcoming.*/
/* Empirical application - macro */
/* Table 5 - main paper */

p_val=0.10;
p_val1=0.05;
p_test=0.05;
delta=1/4;
delta1=1/2;
pmax=4;

a_size=0.10;
a_size1=0.05;
delta_thr=1/2;
delta2_thr=1/4;

all_n={};
all_a14={};
all_acip14={};
all_acim14={};
all_a12={};
all_acip12={};
all_acim12={};
all_a141={};
all_acip141={};
all_acim141={};
all_a121={};
all_acip121={};
all_acim121={};

/* PART I: UPLOAD THE DATA */
/* Load Data */
x_swfq=xlsreadm("FRED_QSWf.xls","b2:gf127",1,0); 

/*x=x_swfq[1:80,.];*/ /* Q1:1988-Q4:2007 */
x=x_swfq; /* Q1:1988-Q2:2019 */

t=rows(x); 
n=cols(x);
    
/* Regression with CSA or observed factors */
  ff=meanc(x');
  {a_hat14p,rej,ci_p14,ci_m14,mR2_14}=a_demp(x,ff,p_val,delta,p_test);
  {a_hat12p,rej,ci_p12,ci_m12,mR2_12}=a_demp(x,ff,p_val,delta1,p_test);
  {a_hat14p1,rej,ci_p141,ci_m141,mR2_141}=a_demp(x,ff,p_val1,delta,p_test);
  {a_hat12p1,rej,ci_p121,ci_m121,mR2_121}=a_demp(x,ff,p_val1,delta1,p_test);
  a_14=a_hat14p[1,cols(a_hat14p)];
  a_cip14=ci_p14[1,cols(a_hat14p)];
  a_cim14=ci_m14[1,cols(a_hat14p)];  
  a_12=a_hat12p[1,cols(a_hat12p)];  
  a_cip12=ci_p12[1,cols(a_hat12p)];
  a_cim12=ci_m12[1,cols(a_hat12p)];
  a_141=a_hat14p1[1,cols(a_hat14p1)];
  a_cip141=ci_p141[1,cols(a_hat14p1)];
  a_cim141=ci_m141[1,cols(a_hat14p1)];  
  a_121=a_hat12p1[1,cols(a_hat12p1)];  
  a_cip121=ci_p121[1,cols(a_hat12p1)];
  a_cim121=ci_m121[1,cols(a_hat12p1)];    
  
/* Exponent of CSD */  
  x_st=standx(x);
  {rej,a_tilde,v_f_tilde,rej,rej,ggg_tilde,s_tilde,rej}=atildeall(x_st,a_size);
  om_tilde=1.65*((1/t)*(v_f_tilde)+(4/n)*(n^(1-a_tilde)*s_tilde/(ggg_tilde-1)))^(1/2)/(2*ln(n)); /* 90% CI: 1.65, 95% CI: 1.96 */
  {rej,a_tilde1,v_f_tilde1,rej,rej,ggg_tilde1,s_tilde1,rej}=atildeall(x_st,a_size1);
  om_tilde1=1.65*((1/t)*(v_f_tilde1)+(4/n)*(n^(1-a_tilde1)*s_tilde1/(ggg_tilde1-1)))^(1/2)/(2*ln(n)); /* 90% CI: 1.65, 95% CI: 1.96 */
  /* Save */
  all_n=all_n|n;
  all_a14=all_a14|a_14;
  all_acip14=all_acip14|a_cip14;
  all_acim14=all_acim14|a_cim14;
  all_a12=all_a12|a_12;
  all_acip12=all_acip12|a_cip12;
  all_acim12=all_acim12|a_cim12;
  all_a141=all_a141|a_141;
  all_acip141=all_acip141|a_cip141;
  all_acim141=all_acim141|a_cim141;
  all_a121=all_a121|a_121;
  all_acip121=all_acip121|a_cip121;
  all_acim121=all_acim121|a_cim121;
  
"n";all_n;
"t";t;
"Strength of CSA (delta=1/4), Strength of CSA (delta=1/2), Exponent of CSD";
"p=0.10";
all_acim14~all_a14~all_acip14;
all_acim12~all_a12~all_acip12;
(a_tilde-om_tilde)~a_tilde~(a_tilde+om_tilde);
"p=0.05";
all_acim141~all_a141~all_acip141;
all_acim121~all_a121~all_acip121;
(a_tilde1-om_tilde1)~a_tilde1~(a_tilde1+om_tilde1);

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


/* a_tilde estimates: temporal structure on factors; cross-sectional dependence on errors */
proc (8)=atildeall(x,a_size);
local n,t,p,z,ln_z,x_bar1,x_bar1_c,std_x_bar1,pc,ln_var,x_bar,m_x_bar,x_bar_stand,i,x_bar_2m,m_x_bar_2m,x_bar_2m_st,x_bar_2m_st_lag,
v_all,v_1,dv,rhs,b,s_b,e_nw,sse_nw,sig2_nw,v_f_2,c_avg,m_c_avg,c_avg_2,index,e,e_bar,m_e_bar,e_bar_stand,e_bar_2m,v_e_bar,s_hat,a_dot,
a_tilde,rhsx,coefx,residx,ssex,sdx,dfx,t_test,p_n,theta,size,x_str,x_str1,xstr_bar,musqr_thr,a_thrtilde,ggg_o,c_avg_sel,m_c_avg_sel,
frasel,s_frasel,ggg_t,c_avg_selt,m_c_avg_selt,fraselt,s_fraselt,j,s_ttest,s_size;

n=cols(x);
t=rows(x);
p=ceil(t^(1/3));
z=seqa(1,1,n); /* cross sectional trend */
ln_z=ln(z);
x_bar1=meanc(x');
x_bar1_c=x_bar1;
std_x_bar1=stdc(x_bar1_c);
ln_var=ln(std_x_bar1^2);
x_bar=x_bar1_c./std_x_bar1; /* standardise the cross-sectional avgs */
m_x_bar=meanc(x_bar);
x_bar_stand=zeros(t,1);
i=1;
do while i<=t;
x_bar_stand[i,.]=x_bar[i,.]-m_x_bar;
i=i+1;
endo;
/*Newey-West method*/
x_bar_2m=x_bar_stand^2;
m_x_bar_2m=meanc(x_bar_2m);
x_bar_2m_st=zeros(t,1);
for i(1,t,1);
x_bar_2m_st[i,.]=x_bar_2m[i,.]-m_x_bar_2m;
endfor;
x_bar_2m_st_lag=zeros(rows(x_bar_2m_st),p);
for i(1,p,1);
x_bar_2m_st_lag[.,i]=lagn(x_bar_2m_st,i);
endfor;
v_all=x_bar_2m_st~x_bar_2m_st_lag;
v_1=v_all[p+1:t,.];
dv=v_1[.,1];
rhs=v_1[.,2:p+1];
b=inv(rhs'*rhs)*rhs'*dv;
s_b=sumc(b);
e_nw=dv-rhs*b;
sse_nw=e_nw'*e_nw;
sig2_nw=sse_nw/(t-cols(rhs));
v_f_2=sig2_nw/(1-s_b)^2;

c_avg=inv(x_bar'*x_bar)*x_bar'*x; /* OLS estimate standardised cross-sectional coefficients*/
m_c_avg=meanc(c_avg');
{pc}=getpc(x,4);
c_avg_2=inv(pc'*pc)*pc'*x; /* OLS estimate non standardised cross-sectional coefficients*/
index=rev(sortc(abs(c_avg_2')~z,1));
e=x-pc*c_avg_2; /* calculate residuals from non standardised cross-sectionals regression */
e_bar=meanc(e');
m_e_bar=meanc(e_bar);
e_bar_stand=zeros(t,1);
i=1;
do while i<=t;
e_bar_stand[i,.]=e_bar[i,.]-m_e_bar;
i=i+1;
endo;
e_bar_2m=e_bar_stand^2;
v_e_bar=meanc(e_bar_2m);
s_hat=n*v_e_bar;

a_dot=1+(1/2)*(ln_var/ln(n));
a_tilde=a_dot-(1/2)*(s_hat/(n*ln(n)*std_x_bar1^2));

rhsx=ones(t,1)~x_bar1;
coefx=inv(rhsx'*rhsx)*rhsx'*x;
residx=x-rhsx*coefx;
ssex=residx'*residx/(t-cols(rhsx));
sdx=sqrt(diag(ssex*inv(x_bar1'*x_bar1)));
dfx=t-cols(rhsx);
t_test=coefx[2,.]'./sdx[.,1];
size=zeros(n,1);
x_str=zeros(t,n);
s_ttest=rev(sortc(abs(t_test)~z,1));
j=1;
do while  j<=cols(coefx);
p_n=a_size/(n-j+1);
theta=cdfni(1-p_n/2);
  if abs(s_ttest[j,1])>=theta;
   size[j,1]=1;
  else;
   size[j,1]=0;
  endif;
j=j+1;
endo;
s_size=sortc(size~s_ttest[.,2],2);
x_str=s_size[.,1]'.*x;
x_str=x_str';
x_str1=delif(x_str,x_str[.,1].==0);
if x_str1==miss(1,1);
   musqr_thr=1;
  else;
   x_str1=x_str1';
   xstr_bar=meanc(x_str1');
   musqr_thr=meanc((xstr_bar-meanc(xstr_bar))^2);
endif;
a_thrtilde=a_tilde-(1/2)*(ln(musqr_thr)/ln(n));

ggg_o=round(n^(a_tilde));
    if ggg_o>=n;
       ggg_o=n;
    elseif ggg_o<1;
        ggg_o=1;
    else;
       ggg_o=ggg_o;
    endif;
c_avg_sel=rev(sortc(c_avg'~abs(c_avg'),2));
c_avg_sel=c_avg_sel[1:ggg_o,1];
m_c_avg_sel=meanc(c_avg_sel);
frasel=zeros(1,ggg_o);
for i(1,ggg_o,1);
  frasel[.,i]=(c_avg_sel[i,1]-m_c_avg_sel)^2;
endfor;
s_frasel=sumc(frasel');

ggg_t=round(n^(a_thrtilde));
    if ggg_t>=n;
       ggg_t=n;
    elseif ggg_t<1;
        ggg_t=1;
    else;
       ggg_t=ggg_t;
    endif;
c_avg_selt=rev(sortc(c_avg'~abs(c_avg'),2));
c_avg_selt=c_avg_selt[1:ggg_t,1];
m_c_avg_selt=meanc(c_avg_selt);
fraselt=zeros(1,ggg_t);
for i(1,ggg_t,1);
  fraselt[.,i]=(c_avg_selt[i,1]-m_c_avg_selt)^2;
endfor;
s_fraselt=sumc(fraselt');
retp(a_tilde,a_thrtilde,v_f_2,ggg_o,s_frasel,ggg_t,s_fraselt,musqr_thr); /*output*/
endp;


proc getpc(ax,j1);
local j,x,n,t,first_j,eigvals,eigvecs,evals,evecs,pc;
x=ax;
n=cols(ax);
t=rows(ax);
j=j1; /* or 'factor' if more than the maximum pc is needed */
if n<t;
first_j=seqa(n,-1,j);
{eigvals,eigvecs}=eigrs2(x'x); /*sorted in ascending order*/
evals=submat(eigvals,first_j,1); /* pick j largest eigenvalues */
evecs=submat(eigvecs,0,first_j);
pc=x*evecs;
elseif n>=t;
first_j=seqa(t,-1,j);
{eigvals,eigvecs}=eigrs2(x*x'); /*sorted in ascending order*/
evals=submat(eigvals,first_j,1); /* pick j largest eigenvalues */
evecs=submat(eigvecs,0,first_j);
pc=evecs; /* eigenvectors of xx' = PC of x */
endif;
retp(pc);
endp;


/* Standardise data */
proc standx(x);
local n,t,m_x,std_x,x_stand;
t=rows(x);
n=cols(x);
m_x=meanc(x)';
std_x=stdc(x)';
x_stand=zeros(t,n);
for i(1,n,1);
x_stand[.,i]=(x[.,i]-m_x[1,i])/std_x[1,i];
endfor;
retp(x_stand);
endp;
