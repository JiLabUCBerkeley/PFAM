function zernikeCoeff = PTT2ZernikeCoeff(PTT)

load('data/ZernikeCoeff2PTT2.mat', 'PTT_mode'); 

ZernikeModeN = 55;
zernikeCoeff = zeros(1, ZernikeModeN);
for i = 1:ZernikeModeN
    Mode = PTT_mode(i, :, :);
    zernikeCoeff(i) = sum(PTT(:) .* Mode(:)) / sum(Mode(:).^2);   % Zernike_coef*PTT_mode=PTT --> Zer_coef=PTT/mode 
end

end

