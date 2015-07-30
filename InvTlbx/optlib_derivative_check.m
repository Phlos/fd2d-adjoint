function [dcheck] = optlib_derivative_check(m,dm,hpmin,hpmax,usr_par)
%
%

[j, g] = eval_objective_and_gradient(m, optlib_generate_random_string(8), usr_par);
j
norm(g)
norm_dm = norm(dm)

dcheck = zeros(hpmax-hpmin+1,5);
it=0;
jh_vec=zeros(hpmax - hpmin + 1, 1);
for hp=hpmin:1:hpmax
    it=it+1;
    mh = m + 10^hp * dm;
    [jh] = eval_objective(mh,optlib_generate_random_string(8), usr_par);
    jh_vec(it) = jh;
    djdm = g' * dm;
    djdmh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdm, djdmh, abs(djdm - djdmh), abs(djdm - djdmh) / abs(djdm)];
end

h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
j
jh_vec

end