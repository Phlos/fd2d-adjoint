function v_out = make_seismogram_zeros(v_in)

% function to make a seismogram struct of the same shape as v_in, but zeros

for ii = 1:length(v_in)
    
    comps = fieldnames(v_in{ii});
    for jj = 1:length(comps)
        v_out{ii}.(comps{jj}) = zeros( size(v_in{ii}.(comps{jj})) );
    end
end

end