%% norm01.m
%% Author: Erin Gilbert
%% Created: May 24 2018
%% Modified: Sept 16 2018

%% Usage: norm01(ARRAY_OF_VALUES)
%% Purpose: Normalize data [0,1]

function [n] = norm01(x)
    n= (x - min(min(x)))./(max(max(x))-min(min(x)));
end
