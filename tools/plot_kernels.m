
function fig_knl = plot_kernels(K)

% plots absolute kernels

%- initialise input -------------------------------------------------------
path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);


%- plot -------------------------------------------------------------------

% figure position and size

fn_params = fieldnames(K);

count = 0;
for i = 1:length(fn_params)
    fn_comps = fieldnames(K.(fn_params{i}));
    for j = 1:length(fn_comps)
        if ((i==1) && (j==1))
            kernels_together = K.(fn_params{i}).(fn_comps{j});
        else
            kernels_together = kernels_together + K.(fn_params{i}).(fn_comps{j});
        end
        count = count+1;
    end
end

scale = prctile(abs(kernels_together(:)),99.0);

fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)
                


% actual plotting
% count = 0;
for i = 1:length(fn_params)
    fn_comps = fieldnames(K.(fn_params{i}));
    for j = 1:length(fn_comps)
%         count = count+1;
        subplot(3,ceil(nrec/3),(i-1)*3 + j)
        whichkernel = [fn_params{i}, ' ', fn_comps{j}];
        kernel = K.(fn_params{i}).(fn_comps{j});
        
        % getting potential +Inf and -Inf out of the system
        minnoninf = min(isfinite(kernel(:)) .* kernel(:) );
        maxnoninf = max(isfinite(kernel(:)) .* kernel(:) );
        kernel(kernel==-Inf) = minnoninf;
        kernel(kernel==Inf) = maxnoninf;
    
        scale = prctile(abs(kernel(:)),99);
        plot_kernel(X,Z,kernel,whichkernel,'fixed',scale,stf_PSV);
%         scale = prctile(abs(kernel(:)),99)
    end
end

% subplot (3,2,1);
% plot_kernel(X,Z,K_rel.rho2.PSV,'relative density PSV','fixed',scale,stf_PSV);
% subplot (3,2,3);
% plot_kernel(X,Z,K_rel.vs2.PSV,'relative v_s PSV','fixed',scale,stf_PSV);
% subplot (3,2,5);
% plot_kernel(X,Z,K_rel.vp2.PSV,'relative v_p PSV','fixed',scale,stf_PSV);
% 
% subplot (3,2,2);
% plot_kernel(X,Z,K_rel.rho2.SH,'relative density SH','fixed',scale,stf_PSV);
% subplot (3,2,4);
% plot_kernel(X,Z,K_rel.vs2.SH,'relative v_s SH','fixed',scale,stf_PSV);
% subplot (3,2,6);
% plot_kernel(X,Z,zeros(size(X')),'relative v_p SH','perc',100,stf_PSV);

caxis([-scale scale]);
h = colorbar('horiz');
set(h,'Position',[0.3 0.05 0.4 0.02])