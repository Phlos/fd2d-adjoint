function [sig,fn,gn]=optlib_stepsize_wolfe(xj,s,stg,f,delta,theta,sig0,usr_par)
%
%

s=-s;
stg=-stg;

sig=sig0;
 xn=xj-sig*s;
 [fn] = eval_objective(xn, usr_par);
 %[fn]=feval(fct,xn);
% Determine maximal sig=sig0/2^k satisfying Armijo
 while (f-fn<delta*sig*stg)
  sig=0.5*sig;
  xn=xj-sig*s;
  [fn] = eval_objective(xn, usr_par);
  %[fn]=feval(fct,xn);
 end
 [gn] = eval_grad_objective(xn, usr_par);

% If sig=sig0 satisfies Armijo then try sig=2^k*sig0
% until sig satisfies also the Wolfe condition
% or until sigp=2^(k+1)*sig0 violates the Armijo condition
 if (sig==sig0)
  xnn=xj-2*sig*s;
  [fnn,gnn] = eval_objective_and_gradient(xnn, usr_par);
  %[fnn,gnn]=feval(fct,xnn);
  while (gn'*s>theta*stg)&(f-fnn>=2*delta*sig*stg)
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
  [fnn,gnn] = eval_objective_and_gradient(xnn, usr_par);
  %[fnn,gnn]=feval(fct,xb);
  if (f-fnn>=delta*sigb*stg)
   sig=sigb;
   xn=xb;
   fn=fnn;
   gn=gnn;
  else
   sigp=sigb;
  end
 end


% % % 
% % %  sig=sig0;
% % %  xn=xj+sig*s;
% % %  [fn] = eval_objective(xn, usr_par);
% % % % Determine maximal sig=sig0/2^k satisfying Armijo
% % %  while (fn-f>delta*sig*stg)
% % %   sig=0.5*sig;
% % %   xn=xj+sig*s;
% % %   [fn] = eval_objective(xn, usr_par);
% % %  end
% % %  [gn] = eval_grad_objective(xn, usr_par);
% % % 
% % % % If sig=sig0 satisfies Armijo then try sig=2^k*sig0
% % % % until sig satisfies also the Wolfe condition
% % % % or until sigp=2^(k+1)*sig0 violates the Armijo condition
% % %  if (sig==sig0)
% % %   xnn=xj+2*sig*s;
% % %   [fnn, gnn] = eval_objective_and_gradient(xnn, usr_par);
% % %   %[fnn,gnn]=feval(fct,xnn);
% % %   while (gn'*s  theta * stg) & (fnn-f<=2*delta*sig*stg)
% % %    sig=2*sig;
% % %    xn=xnn;
% % %    fn=fnn;
% % %    gn=gnn;
% % %    xnn=xj+2*sig*s;
% % %    [fnn, gnn] = eval_objective_and_gradient(xnn, usr_par);
% % %   end
% % %  end
% % %  sigp=2*sig;
% % % 
% % % % Perform bisektion until sig satisfies also the Wolfe condition
% % %  while (gn'*s<theta*stg)
% % %   sigb=0.5*(sig+sigp);
% % %   xb=xj+sigb*s;
% % %   [fnn, gnn] = eval_objective_and_gradient(xb, usr_par);
% % %   %[fnn,gnn]=feval(fct,xb);
% % %   if (fnn - f <= delta*sigb*stg)
% % %    sig=sigb;
% % %    xn=xb;
% % %    fn=fnn;
% % %    gn=gnn;
% % %   else
% % %    sigp=sigb;
% % %   end
 end
