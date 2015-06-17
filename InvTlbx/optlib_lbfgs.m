function [flag,mfinal]=optlib_lbfgs(m0, usr_par, initial_steplength, tolerance, max_iterations, output_file)
%
%

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % constants 0<del<theta<1, del<1/2 for strong Wolfe condition
    delta=0.001;
    theta=0.6;
    % constant 0<al<1 for sufficient decrease condition
    alpha=0.001;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  


fid = fopen(output_file,'a+');
it=0;

model.m = m0;
model.string = optlib_generate_random_string(8);

[model.objective, model.gradient] = eval_objective_and_gradient(model.m, model.string, usr_par);

normg0=norm(model.gradient);
model.normg=normg0;

fprintf(fid,'it=%d   j=%e   ||g||=%e  \n',it,model.objective,model.normg);


% init data for LBFGS here
[previous_models, previous_gradients] = optlib_restore_lbfgs_information(usr_par);

if (size(previous_models,2) > 0)
   l_history = size(previous_models,2);
   lmax = max(5,l_history);
    
    P=zeros(size(m0,1),lmax);
    D=zeros(size(m0,1),lmax);
    ga=zeros(lmax,1);
    rho=zeros(lmax,1);
    
    for i=2:l_history
       P(:,i-1) = previous_models(:,i-1) - previous_models(:,i);
       D(:,i-1) = previous_gradients(:,i-1) - previous_gradients(:,i);
    end
    P(:,l_history) = previous_models(:,l_history) - model.m;
    D(:,l_history) = previous_gradients(:,l_history) - model.gradient;
    
    for i=1:l_history
       rho(i) = 1.0 / (P(:,i)' * D(:,i)); 
    end
    gak = P(:,l_history)' * D(:,l_history) / norm(P(:,l_history));
    l=min(lmax, l_history);
    ln=min(lmax, l_history);

else
    lmax=5;
    P=zeros(size(m0,1),lmax);
    D=zeros(size(m0,1),lmax);
    ga=zeros(lmax,1);
    rho=zeros(lmax,1);
    l=0;
    ln=1;
    gak=1;
end


% init data for LBFGS here

fid = fopen(output_file,'a+');
it=0;

model.m = m0;
model.string = optlib_generate_random_string(8);

[model.objective, model.gradient] = eval_objective_and_gradient(model.m, model.string, usr_par);

normg0=norm(model.gradient);
model.normg=normg0;

fprintf(fid,'it=%d   j=%e   ||g||=%e  \n',it,model.objective,model.normg);

% main loop
while (model.normg>tolerance*normg0 && it < max_iterations)
    
    it=it+1;

    % compute BFGS-step s=B*g;
    g=model.gradient;
    q=g;
    for j=1:l
        i=mod(ln-j-1,lmax)+1;
        ga(i)=rho(i)*(P(:,i)'*q);
        q=q-ga(i)*D(:,i);
    end
    r=gak*q;
    for j=l:-1:1
        i=mod(ln-j-1,lmax)+1;
        be=rho(i)*(D(:,i)'*r);
        r=r+(ga(i)-be)*P(:,i);
    end
    s=r;
    step='LBFGS';

    % check if BFGS-step provides sufficient decrease; else take gradient
    stg=s'*g;
    if stg<min(alpha,model.normg)*model.normg*norm(s)
        s=g;
        stg=s'*g;
        step='Grad';
    end

    % choose sigma by Powell-Wolfe stepsize rule
    if (it==1)
        % TODO: provide alternative guess of the initial step length based
        % on relative model perturbations
        sigma=initial_steplength;
    else
        sigma = 1.0;
    end

    % TODO: use quadratic approximation to compute step-length 
    % for iteration 1
    %[sigma,objective_new,gn]=optlib_wolfe(m,s,stg,objective,delta,theta,sigma,usr_par);
    [sigma,model_new]=optlib_wolfe(model.m,s,stg,model.objective,delta,theta,sigma,usr_par);
 
    mn=model_new.m;
    gn=model_new.gradient;
    model_new.normg = norm(model_new.gradient);
    
    fprintf(fid,'it=%d   j=%e   ||g||=%e   sigma=%6.5f ||s||=%e  step=%s\n', ...
            it, model_new.objective, model_new.normg,sigma,norm(s),step);
        
    % update BFGS-matrix
    d=g-gn;
    p=model.m-mn;
    dtp=d'*p;
    if dtp>=1e-8*norm(d)*norm(p)
        rho(ln)=1/dtp;
        D(:,ln)=d;
        P(:,ln)=p;
        l=min(l+1,lmax);
        ln=mod(ln,lmax)+1;
        %  if l==lmax
            gak=dtp/(d'*d);
        % end
    else
        fprintf(fid, '\tWarning: Bad angle of search direction. Update of BFGS curvature information skipped.\n');
    end
    
    model = model_new;
    
    usr_par = new_iteration(it, model.m, model.objective, model.gradient, usr_par);

end


if (model.normg<=tol*normg0)
    fprintf(fid,'Successful termination with ||g||<%e*min(1,||g0||):\n',tolerance);
    flag = 0;
else
    fprintf(fid,'Maximum number of iterations reached.\n');
    flag = 1;
end

% return final model and flag
mfinal=m;
fclose(fid);
