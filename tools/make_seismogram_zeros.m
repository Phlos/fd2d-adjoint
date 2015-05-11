function v_out = make_seismogram_zeros(v_in)

% function to make a 'seismogram' struct of the same shape as v_in, but zeros
% shape of v_in: v_in{irec}.comp(jj) -- v_in{1}.x/y/z    (cell)
%            or: v_in(irec).comp(jj) -- v_in(1).x/y/z    (struct)

if iscell(v_in)
    for ii = 1:length(v_in)
        
        comps = fieldnames(v_in{ii});
        for jj = 1:length(comps)
            v_out{ii}.(comps{jj}) = zeros( size(v_in{ii}.(comps{jj})) );
        end
    end
elseif isstruct(v_in)
    for ii = 1:length(v_in)
        comps = fieldnames(v_in(ii));
        for jj = 1:length(comps)
            v_out(ii).(comps{jj}) = zeros( size(v_in(ii).(comps{jj})) );
        end
    end
end

end