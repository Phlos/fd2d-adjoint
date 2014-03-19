% plot src and rec in adjoint field

for k=1:ns
    plot(src_x(k),src_z(k),'ko');
end
for k=1:length(orig_src_x)
    plot(orig_src_x(k),orig_src_z(k),'kx');
end