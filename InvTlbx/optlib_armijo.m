function [sig,fn]=optlib_armijo(xj,s,stg,f,delta,sig0,usr_par)
%
%
%
 sig=sig0;
 xn=xj+sig*s;
 [fn] = eval_objective(xn, usr_par);
% determine maximal sig=sig0/2^k satisfying Armijo
 while (fn-f>delta*sig*stg)
  sig=0.5*sig;
  xn=xj+sig*s;
  [fn] = eval_objective(xn, usr_par);
 end

% if sig=sig0 satisfies Armijo and sig_gt_1 is given then try sig=2^k*sig0
 if (sig==sig0) %& (nargin==8)
  xnn=xj+2*sig*s;
  [fnn]= eval_objective(xnn, usr_par);
  while (fnn-f<=2*delta*sig*stg)
   sig=2*sig;
   xn=xnn;
   fn=fnn;
   xnn=xj+2*sig*s;
   fnn=eval_objective(xnn, usr_par);
  end
 end
