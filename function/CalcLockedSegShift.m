function Values = CalcLockedSegShift(Values, Segs, Param, Convert_hex2jk, Convert_jk2hex)

for i =  Segs
    % disp(i)
    [neighbours_hex, weights_hex] = Find6NeighbrIdx(i, Param, Convert_hex2jk, Convert_jk2hex);
    for j = 1:size(Values,2)
        Values(i,j) = mean(Values(neighbours_hex, j) .* weights_hex);
    end
end

end