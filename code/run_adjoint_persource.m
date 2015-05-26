function [Kseis, sEventKnls] = run_adjoint_persource(Model, sEventAdstf)

% wrapper to run adjoint for each kernel consecutively

%% prep
input_parameters;
nsrc = length(sEventAdstf);

%% actual loop over sources

% per source:   
for isrc = 1:nsrc
    disp(['Running adjoint wave propagation - src. nr. ',num2str(isrc),'/',num2str(nsrc)]);
    
    % obtain correct adstf
    adstf = sEventAdstf(isrc).adstf;
    
    % load u_fw, v_fw from file
    prevmsg = sprintf('loading forward fields...');
    fprintf(prevmsg);
    load(['../output/forwardfield.src-',num2str(isrc),'.mat'], 'u_fw', 'v_fw');
    reverseStr = repmat(sprintf('\b'), 1, length(prevmsg));
    fprintf(reverseStr);
    
    % run adjoint with those u_fw, v_fw
    Kevent = run_adjoint(u_fw, v_fw, adstf, Model);
        
    % save adjoint into sEventKnls
%     sEventKnls(isrc) = Kevent;
    sEventKnls(isrc).rho.total = Kevent.rho.total;
    sEventKnls(isrc).mu.total = Kevent.mu.total;
    sEventKnls(isrc).lambda.total = Kevent.lambda.total;
    if strcmp(wave_propagation_type, {'both', 'SH'})
        sEventKnls(isrc).rho.SH = Kevent.rho.SH;
        sEventKnls(isrc).mu.SH = Kevent.mu.SH;
        sEventKnls(isrc).lambda.SH = Kevent.lambda.SH;
    end
    
end

% add up all event kernels into a single big one
Kseis = add_kernels(sEventKnls);



end

function Kadded = add_kernels(Kin)

nsrc = length(Kin);

% initialise Kadded
Ksingle = Kin(1);
comps = fieldnames(Ksingle);
for ii = 1:length(comps)
    fields = fieldnames(Ksingle.(comps{ii}));
    for jj = 1:length(fields);
        Kadded.(comps{ii}).(fields{jj}) = zeros(size(Ksingle.(comps{ii}).(fields{jj})));
    end
end


% fill Kadded
for isrc = 1:nsrc
%     disp(['adding kernel ',num2str(isrc)]);
    Ksingle = Kin(isrc);
    comps = fieldnames(Ksingle);
    for ii = 1:length(comps)
        fields = fieldnames(Ksingle.(comps{ii}));
        for jj = 1:length(fields);
                Kadded.(comps{ii}).(fields{jj}) = ...
                    Kadded.(comps{ii}).(fields{jj}) + Ksingle.(comps{ii}).(fields{jj});

        end
    end

end

end