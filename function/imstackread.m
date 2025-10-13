function ImStack = imstackread(FileName)
ImInfo = imfinfo(FileName);
N = length(ImInfo);
ImStack = uint16(zeros(ImInfo(1).Height, ImInfo(1).Width, length(ImInfo)));
% ImStack = uint16(zeros(ImInfo(1).Width, ImInfo(1).Height, length(ImInfo)));
for ii = 1:N
    ImStack(:,:,ii) = imread(FileName, 'Index', ii, 'Info', ImInfo);
end
