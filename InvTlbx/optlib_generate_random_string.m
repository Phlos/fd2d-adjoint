function [s] = optlib_generate_random_string(n)
% OPTLIB_GENERATE_RANDOM_STRING function to generate random string of
% length n. The string contains letters and numbers.
%
% INPUT:
% n : length of random string
%
% OUTPUT:
% s : random string of length n

    symbols = ['a':'z' 'A':'Z' '0':'9'];
    nums = randi(numel(symbols),[1 n]);
    s = symbols (nums);
end
