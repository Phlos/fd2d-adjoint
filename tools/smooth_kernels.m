function Ksmooth = smooth_kernels(kernel, npoints)

% this function smooths the sensitivity kernels so that no odd
% singularities are formed -- what normally happens close to sources and
% receivers.

%- initialising
path(path,'../input');
input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

disp 'Smoothing kernels!!'

%- gaussian smoothing filter
filt = fspecial('gaussian',[npoints npoints],5);   %  (fspecial is normally part of the image processing toolbox but I adapted
                                                   %  the (probably exactly equivalent) Octave code. fspecial found in tools/)
                                                   % not sure what the 2 means

%-- smooth all 'total' kernels (in rho mu lambda parametrisation)
bips = figure;
for sname = {'rho' 'mu' 'lambda'}
    Ksmooth.(sname{1}) = conv2(kernel.(sname{1}).total, filt, 'same');
    
    % plot the smoothed kernels
    figure(bips);
    cmax = prctile(abs(kernel.(sname{1}).total(:)),99.97);
    subplot(2,1,1)
    knl_title = ['original ',sname{1}];
    plot_kernel(X,Z,kernel.(sname{1}).total,knl_title,'fixed',cmax,[1 0]);
    colorbar;
    subplot(2,1,2)
    knl_title = ['smoothed ',sname{1}];
    plot_kernel(X,Z,Ksmooth.(sname{1}),knl_title,'fixed',cmax,[1 0]);
    colorbar;
%     pause(3.0);
    pause(0.5);
    
end

close(clf);

end