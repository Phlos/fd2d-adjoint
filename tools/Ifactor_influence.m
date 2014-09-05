function Ifactor_influence(dm_size);
% testing at which distance a mass dm has the biggest influence on the ifactor


Mrod = 1;
Lrod = 1;
Irod = 1/3 * Mrod * Lrod^2;
dm = dm_size*Mrod;
r = 0:0.01:1;
Idm = dm*r.^2;
Itot = Irod + Idm;
Ifactor_rod = Irod / (Mrod*Lrod^2);
Ifactor_tot = Itot ./ ( (Mrod + dm) * Lrod^2);
figure
plot(r,Ifactor_tot);
Irod_array = Irod * ones(1,size(r,2));
hold on
plot(r,Irod_array)
