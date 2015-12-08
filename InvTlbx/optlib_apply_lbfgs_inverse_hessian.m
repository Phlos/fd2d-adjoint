
function [BV] = optlib_apply_lbfgs_inverse_hessian(V, LBFGS_data)

    nvec = size(V,2);
    %ga = zeros(LBFGS_data.lmax,1);

    BV = zeros(size(V,1),nvec);
% for nv=1:nvec
%     q=V(:,nv);
%     for j=1:LBFGS_data.l
%         i=mod(LBFGS_data.ln-j-1,LBFGS_data.lmax)+1;
%         ga(i)=LBFGS_data.rho(i)*(LBFGS_data.P(:,i)'*q);
%         q=q-ga(i)*LBFGS_data.D(:,i);
%     end
%     r=gak*q;
%     for j=LBFGS_data.l:-1:1
%         i=mod(LBFGS_data.ln-j-1,LBFGS_data.lmax)+1;
%         be=LBFGS_data.rho(i)*(LBFGS_data.D(:,i)'*r);
%         r=r+(ga(i)-be)*LBFGS_data.P(:,i);
%     end
%     BV(:,nv)=r;
%     
% end
% BV2=BV;
% 
    ga = zeros(LBFGS_data.lmax,nvec);
    Q=V;
    for j=1:LBFGS_data.l
        i=mod(LBFGS_data.ln-j-1,LBFGS_data.lmax)+1;
        ga(i,:) = LBFGS_data.rho(i) * (LBFGS_data.P(:,i)' * Q);
        Q = Q - LBFGS_data.D(:,i) * ga(i,:);
    end
    BV = LBFGS_data.gak * Q;
    for j=LBFGS_data.l:-1:1
        i=mod(LBFGS_data.ln-j-1,LBFGS_data.lmax)+1;
        be=LBFGS_data.rho(i)*(LBFGS_data.D(:,i)'*BV);
        BV=BV+(ga(i,:)-be)*LBFGS_data.P(:,i);
    end
    
end