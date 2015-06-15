function [sig,fn,gn]=optlib_wolfe(xj,s,stg,f,del,theta,sig0,usr_par)
%
%
% S. Ulbrich, May 2002
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sig,xn,fn]=wolfe(xj,s,stg,fct,f,del,sig0,theta)
%
% Determines stepsize satisfying the Powell-Wolfe conditions
%
% Input:  xj       current point
%         s        search direction (xn=xj-sig*s)
%         stg      stg=s'*g
%         fct      name of a matlab-function [f]=fct(x)
%                  that returns the value of the objective function
%         f        current objective function value f=fct(xj)
%         del      constant 0<del<1/2 in Armijo condition f-fn>=sig*del*stg
%         theta    constant del<theta<1 in Wolfe condition gn'*s<=theta*stg
%         sig0     initial stepsize (usually sig0=1)
%
% Output: sig      stepsize sig satisfying the Armijo condition
%         xn       new point xn=xj-sig*s
%         fn       fn=f(xn)
%
 sig=sig0;
 xn=xj-sig*s;
 [fn] = eval_objective(xn, usr_par);
% Determine maximal sig=sig0/2^k satisfying Armijo
 while (f-fn<del*sig*stg)
  sig=0.5*sig;
  xn=xj-sig*s;
  [fn] = eval_objective(xn, usr_par);
 end
 [gn] = eval_grad_objective(xn, usr_par);

% If sig=sig0 satisfies Armijo then try sig=2^k*sig0
% until sig satisfies also the Wolfe condition
% or until sigp=2^(k+1)*sig0 violates the Armijo condition
 if (sig==sig0)
  xnn=xj-2*sig*s;
  [fnn,gnn] = eval_objective_and_gradient(xnn, usr_par);
  %[fnn,gnn]=feval(fct,xnn);
  while (gn'*s>theta*stg)&(f-fnn>=2*del*sig*stg)
   sig=2*sig;
   xn=xnn;
   fn=fnn;
   gn=gnn;
   xnn=xj-2*sig*s;
   [fnn,gnn] = eval_objective_and_gradient(xnn, usr_par);
   %[fnn,gnn]=feval(fct,xnn);
  end
 end
 sigp=2*sig;

% Perform bisektion until sig satisfies also the Wolfe condition
 while (gn'*s>theta*stg)
  sigb=0.5*(sig+sigp);
  xb=xj-sigb*s;
  [fnn,gnn] = eval_objective_and_gradient(xb, usr_par);
  %[fnn,gnn]=feval(fct,xb);
  if (f-fnn>=del*sigb*stg)
   sig=sigb;
   xn=xb;
   fn=fnn;
   gn=gnn;
  else
   sigp=sigb;
  end
 end
