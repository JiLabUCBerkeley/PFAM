%% three types of index order:
%% 1. hex-index: sequential, same indexing as with the DM segments
%% 2. sq-index: sequential, from left to right, top to bottom
%% 3. jk-index: same order as sq-index, but in Cell structure: 
%%              j is #col, j=-L:L; k is #row in col[j], k(j)=1,...,2*L+1-abs(j)

%% for zonal reconstruction, we use hex-index for most cases because 
%% it simplifies the process of applying corrections to the DM.
%% For phase reconstruction, we use j,k-index to make neighbours location easier.
%% sq-indexing is used as an intermediate index during the convertion, 
%% because sq2hex can be easily read from DM's diagram.

L = 7;
SegN = L*(L+1)*3+1;

%% jk-index --> sq-index
cnt = 0;
IndexCell_jk2sq = cell(2*L+1,1);
for j = -L : L
    Col_j2sq = zeros(2*L+1-abs(j), 1);
    for k = 1 : 2*L+1-abs(j)
        cnt = cnt + 1;
        Col_j2sq(k) = cnt;
    end
    IndexCell_jk2sq{j+L+1} = Col_j2sq;
end

%% sq-index --> jk-index
IndexSeq_sq2jk = zeros(SegN, 2);  % element: (j,k)
cnt = 0;
for j = -L : L
    for k = 1 : 2*L+1-abs(j)
        cnt = cnt + 1;
        IndexSeq_sq2jk(cnt, :) = [j, k];
    end
end


%% sq-index --> hex-indexs
IndexSeq_sq2hex = [ 135:142, ...
                  134, 98:104, 143, ...
                  133, 97, 67:72, 105, 144, ...
                  132, 96, 66, 42:46, 73, 106, 145, ...
                  131, 95, 65, 41, 23:26, 47, 74, 107, 146, ...
                  130, 94, 64, 40, 22, 10:12, 27, 48, 75, 108, 147, ...
                  129, 93, 63, 39, 21, 9, 3:4, 13, 28, 49, 76, 109, 148, ...
                  128, 92, 62, 38, 20, 8, 2, 1, 5, 14, 29, 50, 77, 110, 149, ...
                  169, 127, 91, 61, 37, 19, 7:-1:6, 15, 30, 51, 78, 111, 150, ...
                  168, 126, 90, 60, 36, 18:-1:16, 31, 52, 79, 112, 151, ...
                  167, 125, 89, 59, 35:-1:32, 53, 80, 113, 152, ...
                  166, 124, 88, 58:-1:54, 81, 114, 153, ...
                  165, 123, 87:-1:82, 115, 154, ...
                  164, 122:-1:116, 155, ...
                  163:-1:156 ]';


%%  jk-index --> hex-index                                   
cnt = 0;
IndexCell_jk2hex = cell(2*L+1,1);
for j = -L : L
    Col_j2hex = zeros(2*L+1-abs(j), 1);
    for k = 1 : 2*L+1-abs(j)
        cnt = cnt + 1;
        Col_j2sq = IndexCell_jk2sq{j+L+1};
        Col_j2hex(k) = IndexSeq_sq2hex(Col_j2sq(k));
    end
    IndexCell_jk2hex{j+L+1} = Col_j2hex;
end
% disp(IndexCell_jk2n{L+1});        

%% hex-index --> jk-index
IndexSeq_hex2jk = zeros(SegN, 2);  % element: (j,k)
IndexSeq_hex2sq = zeros(SegN, 1);  
for n_hex = 1:SegN
    n_sq = find(IndexSeq_sq2hex == n_hex);   % hex2sq
    IndexSeq_hex2sq(n_hex) = n_sq;
    IndexSeq_hex2jk(n_hex, :) = IndexSeq_sq2jk(n_sq, :);
end


%%
save('Data\IndexTransfer.mat', 'IndexCell_jk2sq', 'IndexSeq_sq2jk', ...
    'IndexSeq_sq2hex', 'IndexSeq_hex2sq', 'IndexSeq_hex2jk', 'IndexCell_jk2hex');
