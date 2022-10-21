function [STM,Gamma,constant,TR] = extractfromtext(filename1, filename2)

% FUNCTION NAME:
%   extractfromtext
%
% DESCRIPTION:
%   Extracts STM and Gamma from the cpp outcome
%
% INPUT:
%   filename - (string) Name of the data file
%
% OUTPUT:
%   STM - (double []) STM matricies
%   GAMMA - (double []) Gamma matricies 
%   constant - (double []) Constant part of the state
%  eta- double [] Nonlinearity index.

M = readmatrix(filename1);

STM = M(:,1:6); 
Gamma = M(:,7:9);

constant = M(:,10);

TR = readmatrix(filename2);

%TR(4:end) = 0.1*TR(4:end);
end
