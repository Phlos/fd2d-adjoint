
% script to test the kernel magnitude error 
% -- Christian's suggestion with directional derivative --

%% some test properties
% number of test dms i will try
ntest = 5;
% anomaly size (1 = 100% of model)
rel_perturbation = 0.01;
% dm in a random direction or just one point?
dmtype = 'single_point'; % 'random' or 'single_point'

%- standard tee step - may have to be adapted
tee = 1*10 .^ [-1 : -2: -7]; % diverges after 10^-10 --> numerical error?
% tee = 1*10 .^ [-4];

%% prepare starting model, true model and seismic kernel

% %- take a 'real' model -- run forward
% disp 'preparing obs'
% [Obs.model, Obs.v, Obs.t, Obs.props, Obs.g] = prepare_obs(15);
% close all;
%     
% %- take a 'starting' model m -- run forward
% % [Start.model, Start.v, Start.t, Start.props, Start.g] = prepare_obs(10);
% disp '---starting model; run foward'
% Start.model = update_model(10);
% [Start.v,Start.t,u_fw,v_fw,rec_x,rec_z]=run_forward(Start.model); close all;
% t = Start.t;
% 
% %- calculate the misfit between real & starting model ( J(m) )
% %  (can be done for seis and grav both)
% [Start.adstf, Start.misfit_seis] = calc_misfitseis_adstf('waveform_difference',t,Start.v,Obs.v); close all;
% 
% %- calculate kernel(s) gradJ
% disp '---run adjoint'
% normKseis = 1;
% Start.Kseis = run_adjoint(u_fw,v_fw,Start.adstf,'waveform_difference',Start.model,normKseis); close all;

%% prepare directions in which the derivatives are calculated
%- make a 'random' set of dm (in the order of magnitude of realistic 
%  actual variations in the model parameters)
%  randi(2,...) gives pseudorandom integers between 1 and 2.
disp '---determining dee.(parameter) for all'
for ii = 1:ntest
    if strcmp(dmtype,'random')
        for param = fieldnames(Start.model)'
            if strcmp(param{1},'lambda')
                [bincounts, centre] = hist(Start.model.(param{1})(:),100); [~,ix]=max(bincounts); modus.(param{1}) = centre(ix);
                dee{ii}.(param{1}) = ( randi(2,size(Start.model.(param{1})))-1 ) * rel_perturbation * modus.(param{1});
            else
                dee{ii}.(param{1}) = zeros(size(Start.model.(param{1})));
            end
        end
    elseif strcmp(dmtype,'single_point')
        for param = fieldnames(Start.model)'
            dee{ii}.(param{1}) = zeros(size(Start.model.(param{1})));
            if strcmp(param{1},'lambda')
                [bincounts, centre] = hist(Start.model.(param{1})(:),100); [~,ix]=max(bincounts); modus.(param{1}) = centre(ix);
                point = randi( length(dee{ii}.(param{1})(:)) )
                dee{ii}.(param{1})(point) = rel_perturbation * modus.(param{1});
            else
                
            end
        end
    end
        
end

%% for each direction 1:ntest, the dirderiv is calculated using the kernel (L) and the limit of misfits (R).
%  the limit of misfits is calculated by a forward calculation of m + t*dm
%  where t is made progressively smaller (e.g. 10^-4 -- 10^-10)

%- start (outer) loop over different dm
for ii = 1:ntest;
    disp([' --testing perturbed model nr ',num2str(ii)]);
    %-- Left: calculate gradJ^T dm (directional derivative) and save
    dirderiv{ii} = 0;
    for param = fieldnames(dee{ii})'
        dirderiv{ii} = dirderiv{ii} + Start.Kseis.(param{1}).total(:)' * dee{ii}.(param{1})(:);
    end
    disp(['   L: directional derivative = ',num2str(dirderiv{ii})]);
    
    if ~(dirderiv{ii} == 0)
        %-- start (inner) loop over different tee
        %     disp ' --starting loop over tee'
        for jj = 1: length(tee)
            disp(['  -testing step ',num2str(jj),' out of ',num2str(length(tee))]);
            disp(['   (tee = ',num2str(tee(jj)),')']);
            
            %--- calculate perturbed model m + t*dm -- run forward
            pert.model.rho = Start.model.rho + tee(jj)*dee{ii}.rho;
            pert.model.mu = Start.model.mu+ tee(jj)*dee{ii}.mu;
            pert.model.lambda = Start.model.lambda+ tee(jj)*dee{ii}.lambda;
            
            [pert.v,~,~,~,~,~]=run_forward(pert.model);
            
            %--- calculate the misfit between the real & perturbed model ( J(m+t*dm) )
            [~, pert.misfit_seis] = calc_misfitseis_adstf('waveform_difference',t,pert.v,Obs.v);
            
            %--- Right: calculate the limit [J(m+t*dm) - J(m)] / t (could also
            %    be done at the end of the loop)
            limiet{ii}(jj) = (pert.misfit_seis.total - Start.misfit_seis.total) / tee(jj);
            disp(['   R: limiet bij deze tt: ',num2str(limiet{ii}(jj))]);
            disp( '                                 ');
            
            %--- ?devise some test for when the optimal tee has been reached?
            
            %-- end (inner) loop over tee
        end
        
        %-- plot t versus limit
        temporary = ones(length(tee))*dirderiv{ii};
        figure;
        loglog(tee,limiet{ii},tee,temporary);
        semilogx(tee,limiet{ii},tee,temporary);
        xlabel('limit step size')
        ylabel('directional derivative and limit approx')
        title(['test ',num2str(ii)])
        
    end
    
%- end (outer) loop over dm
end 