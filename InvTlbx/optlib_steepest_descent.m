function [flag,mfinal]=optlib_steepest_descent(m0,stepsize, sigma_init ,tol,maxiter, usr_par)
%
%

% constant 0<del<1/2 for Armijo condition
delta=0.001;
% constant del<theta<1 for Wolfe condition
theta=0.6;

sigma=sigma_init;
m=m0;
[j, g] = eval_objective_and_gradient(m, usr_par);

normg0=norm(g);
normg=normg0;

it=0;
fid = fopen('iteration.tab','a+');

fprintf(fid,'it=%d   j=%e   ||g||=%e  \n',it,j,normg);

% main loop
while (normg>tol*normg0 && it < maxiter)
 
 it=it+1;
 sigma0=sigma;
 s=-g;
 stg=s'*g;

 if (stepsize==0)
% choose sigma by Armijo stepsize rule starting with previous
% stepsize sig0. If sig0 is acceptable try sig=2^k*sig0.
  [sigma, jn]= optlib_armijo(m,s,stg,j,delta,sigma0,usr_par);
 else
% choose sigma by Powell-Wolfe stepsize rule starting with previous
% stepsize sig0.
  [sigma, jn, gn]=optlib_stepsize_wolfe(m,s,stg,j,delta,theta,sigma0,usr_par);
 end
 m=m+sigma*s;
 j=jn;
 if (stepsize==0)
     [gn] = eval_grad_objective(m, usr_par);
 end
 g=gn;
 normg=norm(g);
 
 usr_par = new_iteration(m,it,usr_par);
 fprintf(fid,'it=%3.d   f=%e   ||g||=%e   sig=%5.3f\n',it,j,normg,sigma); 
end

if (normg<=tol*normg0)
    fprintf(fid,'Successful termination with ||g||<%e*||g0||:\n',tol);
    flag = 0;
else
    fprintf(fid,'Maximum number of iterations reached.\n');
    flag =1;
end
mfinal=m;
fclose(fid);

end