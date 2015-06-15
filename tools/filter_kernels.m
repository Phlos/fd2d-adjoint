function Kf = filter_kernels(K, parametrisation, sigma)

% filters a kernel of the shape K.[param].total with a Gaussian filter,
% using a smoothing width of sigma grid points in the given parametrisation


switch parametrisation
    case 'rhomulambda'
        Kf.rho.total = filter_2Dfield(K.rho.total, sigma);
        Kf.mu.total = filter_2Dfield(K.mu.total, sigma);
        Kf.lambda.total = filter_2Dfield(K.lambda.total, sigma);
    case 'rhovsvp'
        Kf.rho2.total = filter_2Dfield(K.rho2.total, sigma);
        Kf.vs2.total = filter_2Dfield(K.vs2.total, sigma);
        Kf.vp2.total = filter_2Dfield(K.vp2.total, sigma);
end

end