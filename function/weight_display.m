function weight_dis = weight_display(SH_im, PosInt, Weight, Period, Outliers, Slopes)
if nargin < 6 
    Slopes = zeros(size(PosInt));
end
if nargin < 5
    Outliers = [];
end
    
weight_dis = zeros(size(SH_im));
half_Period = ceil(Period / 2);
Seg_active = 1:size(PosInt,1);
Seg_active(Outliers) = [];
for i = Seg_active
    % hexagon shape
    for y = -half_Period-1 : half_Period + 1
        for x = -(half_Period*2 -abs(y))/sqrt(3)-1 : (half_Period*2-abs(y))/sqrt(3)+1
            weight_dis(round(PosInt(i,2)+y), round(PosInt(i,1)+x)) = ...
                Weight(i) + Slopes(i,1) * x + Slopes(i,2) * y;
        end
    end
end
end