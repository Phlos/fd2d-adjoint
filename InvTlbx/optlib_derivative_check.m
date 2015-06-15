function [dcheck] = optlib_derivative_check(m,dm,hpmin,hpmax,usr_par)
%
%

[j, g] = eval_objective_and_gradient(m, usr_par);
j
norm(g)
norm_dm = norm(dm)

dcheck = zeros(hpmax-hpmin+1,4);
it=0;
for hp=hpmin:1:hpmax
    it=it+1;
    mh = m + 10^hp * dm;
    [jh] = eval_objective(mh, usr_par);
    djdm = g' * dm;
    djdmh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdm, djdmh, abs(djdm - djdmh)];
end

h1 = figure;
loglog(dcheck(:,1),dcheck(:,4))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')

end