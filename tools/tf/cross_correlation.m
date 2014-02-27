%- compute cross correlation function -------------------------------------
%
% function [cc,t_cc]=cross_correlation(f,g,t)
%
% cc_i = sum_j f^*(j) g(j+i)

function [cc,t_cc]=cross_correlation(f,g,t)

%- initialisations --------------------------------------------------------

n=length(f);
cc=zeros(2*n+1);

dt=t(2)-t(1);
t_cc=-n*dt:dt:n*dt;

%- compute brute-force correlation function -------------------------------

for i=-n:n
    for j=1:length(f)
        cc(n+i+1)=cc(n+i+1)+conj(f(j))*g(j+i);
    end
end