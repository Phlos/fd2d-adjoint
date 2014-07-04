function Msmooth = smooth_model(Model, npoints, gwid)

% this function smooths the model parameters Model.{parameter} so that no
% singularities are formed -- what normally happens close to sources and
% receivers in a seismic inversion update.
% Essentially a low-pass filter. Theory in http://www.dspguide.com/ch24/1.htm
% The filter is Gaussian with a width of gwid.
%
% INPUT:
% - Model:  struct with Model.rho, Model.mu and Model.lambda
%           which are nx by nz arrays of the value of the parameters at
%           each point in space inside the modelling domain.
% - npoints:the size of the filter
% - gwid:   width of the gaussian in the filter
%
% OUTPUT:  
% - Msmooth:the smoothed model.


%- initialising
path(path,'../input');
input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

% disp 'Smoothing the model!!'

%- gaussian smoothing filter
filt = fspecial('gaussian', [npoints npoints], gwid);   %  (fspecial is normally part of the image processing toolbox but I adapted
                                                   %  the (probably exactly equivalent) Octave code. fspecial found in tools/)
                                                   % gwid is the width of
                                                   % the gaussian. (could
                                                   % make this input too..)

%-- actual smoothing of the models (in rho mu lambda parametrisation)

% bips = figure;
for sname = {'rho' 'mu' 'lambda'}
    Msmooth.(sname{1}) = conv2(Model.(sname{1}), filt, 'valid');
    
%     % plot the smoothed models?
%     figure(bips);
%     cmax = prctile(abs(Model.(sname{1})(:)),99.97);
%     subplot(2,1,1)
%     knl_title = ['original ',sname{1}];
%     plot_kernel(X,Z,Model.(sname{1}),knl_title,'fixed',cmax,[1 0]);
%     colorbar;
%     subplot(2,1,2)
%     knl_title = ['smoothed ',sname{1}];
%     plot_kernel(X,Z,Msmooth.(sname{1}),knl_title,'fixed',cmax,[1 0]);
%     colorbar;
% %     pause(3.0);
%     pause(0.5);
    
end

% close(clf);

end