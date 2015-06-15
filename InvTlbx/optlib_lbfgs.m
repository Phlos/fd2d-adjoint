function [flag,mfinal]=optlib_lbfgs(m0,tol,maxiter,usr_par)
%
%

% constants 0<del<theta<1, del<1/2 for Wolfe condition
delta=0.001;
theta=0.6;
% constant 0<al<1 for sufficient decrease condition
al=0.001;
lmax=5;
P=zeros(size(m0,1),lmax);
D=zeros(size(m0,1),lmax);
ga=zeros(lmax,1);
rho=zeros(lmax,1);
l=0;
gak=1;

m=m0;
[f, g] = eval_objective_and_gradient(m, usr_par);

normg0=norm(g);
normg=normg0;
it=0;
l=0;
ln=1;

fprintf(1,'it=%d   j=%e   ||g||=%e  \n',it,f,normg);

% main loop
while (norm(g)>tol*max(1,normg0))
 it=it+1;
 sig=1;
% compute BFGS-step s=B*g;
 q=g;
 for j=1:l
  i=mod(ln-j-1,lmax)+1;
  ga(i)=rho(i)*(P(:,i)'*q);
  q=q-ga(i)*D(:,i);
%  pause
 end
 r=gak*q;
 for j=l:-1:1
  i=mod(ln-j-1,lmax)+1;
  be=rho(i)*(D(:,i)'*r);
  r=r+(ga(i)-be)*P(:,i);
 % pause;
 end
 s=r;
 step='LBFGS';

% check if BFGS-step provides sufficient decrease; else take gradient
 stg=s'*g;
 if stg<min(al,normg)*normg*norm(s)
  s=g;
  stg=s'*g;
  step='Grad';
 end
% choose sig by Powell-Wolfe stepsize rule
if (it==1)
sig=1e7;
else
sig = 1.0;
end

 [sig,fn,gn]=optlib_wolfe(m,s,stg,f,delta,theta,sig,usr_par);
 
 mn=m-sig*s;
%plot([m(1),mn(1)],[m(2),mn(2)],'o-')
 fprintf(1,'it=%d   j=%e   ||g||=%e   sig=%6.5f   step=%s\n',it,fn,norm(gn),sig,step);
 usr_par = new_iteration(mn,it,usr_par);
%[fn,gn]=feval(fg,mn);
% update BFGS-matrix
 d=g-gn;
 p=m-mn;
 dtp=d'*p;
 if dtp>=1e-8*norm(d)*norm(p)
  rho(ln)=1/dtp;
  D(:,ln)=d;
  P(:,ln)=p;
  l=min(l+1,lmax);
  ln=mod(ln,lmax)+1;
  if l==lmax
   gak=dtp/(d'*d);
  end
 end
 m=mn;
 g=gn;
 f=fn;
 normg=norm(g);

end


if (normg<=tol*normg0)
    fprintf(1,'Successful termination with ||g||<%e*min(1,||g0||):\n',tol);
    flag = 0;
else
    fprintf(1,'Maximum number of iterations reached.\n');
    flag =1;
end
mfinal=m;
