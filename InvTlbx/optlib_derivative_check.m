function [dcheck, dcheck_struct] = optlib_derivative_check(m,dm,hpmin,hpmax,usr_par)
    
% check if the calculation of a derivative is working correctly
%

[j, g] = eval_objective_and_gradient(m, optlib_generate_random_string(8), usr_par);

j
normg = norm(g)
norm_dm = norm(dm)
djdm = g' * dm


dcheck = zeros(hpmax-hpmin+1,6);
it=0;
jh_vec=zeros(hpmax - hpmin + 1, 1);
for hp=hpmin:1:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, hpmax-hpmin+1);
    mh = m + 10^hp * dm;
    [jh] = eval_objective(mh,optlib_generate_random_string(8), usr_par);
    jh_vec(it) = jh;
    djdm = g' * dm;
    djdmh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdm, djdmh, abs(djdm - djdmh), abs(djdm - djdmh) / abs(djdm), djdm / djdmh];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).djdm_LHS = djdm;
    dcheck_struct(it).djdmh_RHS = djdmh;
    dcheck_struct(it).absdif_LR  = abs(djdm - djdmh);
    dcheck_struct(it).reldif_LR  = abs(djdm - djdmh) / abs(djdm);
    dcheck_struct(it).ratio_LR   = djdm / djdmh;
    
    
end

h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
j
jh_vec

end