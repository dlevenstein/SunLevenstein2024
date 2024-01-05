function [Ca] = CaFromA(A,R,Ca_0,Ca_PSP,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addParameter(p,'linlog','log');

parse(p,varargin{:})

linlog = p.Results.linlog;

%A = (kf_0 + k_CamK.*m)./(kf_0 + k_CamK.*m + kd_0 + k_CaN.*n);
switch linlog
    case 'log'
        Ca = Ca_0 + R.*Ca_PSP.*A;
    case 'lin'
        Ca = Ca_0 + log10(R.*Ca_PSP.*A);
end

        
end

