function [ s ] = Sigmoid( x,x0,k,xmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('x0','var'); x0 = 0; end
if ~exist('k','var'); k = 1; end
if ~exist('xmax','var'); xmax = 1; end

s = xmax./(1+exp(-k.*(x-x0)));


end

