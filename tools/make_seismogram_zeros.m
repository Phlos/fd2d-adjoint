function v_out = make_seismogram_zeros(v_in)

% function to make a seismogram struct of the same shape as v_in, but zeros

for ii = 1:length(v_in)
    
    for comps = fieldnames(v_in{ii})
        v_out{ii}.(comps{1}) = zeros( size(v_in{ii}.(comps{1})) );
    end
end

end