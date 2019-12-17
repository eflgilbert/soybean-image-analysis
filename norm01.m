%Normalize data [0,1]
function [n] = norm01(x)
    n= (x - min(min(x)))./(max(max(x))-min(min(x)));
end
