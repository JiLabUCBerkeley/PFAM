function [neighbours_hex, weights_hex] = Find6NeighbrIdx(SegId, Param, Convert_hex2jk, Convert_jk2hex)

    if nargin < 2
        load('Data\IndexTransfer.mat', 'IndexSeq_hex2jk', 'IndexCell_jk2hex');       % index convertion
        Convert_hex2jk = IndexSeq_hex2jk;
        Convert_jk2hex = IndexCell_jk2hex;
    end

    % determine 6 neighbours surrounding index p (hex-index)        
    neighbours_hex = ones(6,1);    % hex-index, initialize to center segment
    weights_hex = ones(6,1);    % weights of neighbours: 1 if within the pupil
    j = Convert_hex2jk(SegId,1);  % convert hex-index to jk-index
    % disp(SegId)
    % disp(1)

    k = Convert_hex2jk(SegId,2);  % convert hex-index to jk-index
    L = Param.level;
    j = j + L + 1;  % offset j=-L:L to j=1:2L+1
    if j < L+1      % boundary constraints
        if k-1 < 1
            weights_hex(1) = 0;
        else
            neighbours_hex(1) = Convert_jk2hex{j}(k-1);
        end
        if j-1 < 1 || k-1 < 1
                weights_hex(2) = 0;
        else
                neighbours_hex(2) = Convert_jk2hex{j-1}(k-1);
        end
        if j-1 < 1 || k > 2*L+1-abs(j-2-L)
            weights_hex(3) = 0;
        else
            neighbours_hex(3) = Convert_jk2hex{j-1}(k);
        end
        if k+1 > 2*L+1-abs(j-L-1)
            weights_hex(4) = 0;
        else
            neighbours_hex(4) = Convert_jk2hex{j}(k+1);
        end
        neighbours_hex(5) = Convert_jk2hex{j+1}(k+1); % 5,6 are always safe 
        neighbours_hex(6) = Convert_jk2hex{j+1}(k);  

    elseif j == L+1
        if k-1 < 1
            weights_hex(1) = 0;
            weights_hex(2) = 0;
            weights_hex(6) = 0;
        else
            neighbours_hex(1) = Convert_jk2hex{j}(k-1);
            neighbours_hex(2) = Convert_jk2hex{j-1}(k-1);
            neighbours_hex(6) = Convert_jk2hex{j+1}(k-1);
        end
        if k+1 > 2*L+1-abs(j-L-2)
            weights_hex(3) = 0;
        else
            neighbours_hex(3) = Convert_jk2hex{j-1}(k); 
        end
        if k+1 > 2*L+1-abs(j-L-1)
            weights_hex(4) = 0;
        else
            neighbours_hex(4) = Convert_jk2hex{j}(k+1);
        end
        if k+1 > 2*L+1-abs(j-L)
            weights_hex(5) = 0;
        else
            neighbours_hex(5) = Convert_jk2hex{j+1}(k);
        end

    else %j>0
        if k-1 < 1
            weights_hex(1) = 0;
        else
            neighbours_hex(1) = Convert_jk2hex{j}(k-1);
        end
        if j+1 > L*2+1 || k-1 < 1
            weights_hex(6) = 0;
        else
            neighbours_hex(6) = Convert_jk2hex{j+1}(k-1);
        end
        if k+1 > 2*L+1-abs(j-L-1)
            weights_hex(4) = 0;
        else
            neighbours_hex(4) = Convert_jk2hex{j}(k+1);
        end
        if j+1 > L*2+1 || k > 2*L+1-abs(j-L)
            weights_hex(5) = 0;
        else
            neighbours_hex(5) = Convert_jk2hex{j+1}(k);
        end
        neighbours_hex(2) = Convert_jk2hex{j-1}(k);  % 2, 3 are always safe
        neighbours_hex(3) = Convert_jk2hex{j-1}(k+1);
    end

    % set weights of outliers to zero
    for n = 1:6
        if any(Param.Outlier_hex(:) == neighbours_hex(n))
            weights_hex(n) = 0;
        end
    end
end