function [U,D] = IncrementalSKL(A,Uin,Din,K)
%INCREMENTALSKL Function for computing the Sequential Karhunen-Loeve basis
%extraction using the incremental algorithm developed by Levy et Lindenbaum
% En IEEE Trans on Image Proc., 9(8), 2000
%
% [U,S]= IncrementalSKL(A,Uin,Din) computes the left-side U and diagonal D
% matrices from the prexisting Uin, Din matrices when adding a new column
% of data A
%
% @ 2021 Jose A. Gomez-Pedrero, Julio C. Estrada, Juan A. Quiroga, 
% Jose Alonso and Javier Vargas
%
% DISCLAIMER: This software is supplied without any guarantee or 
% any support whatsoever.

if(nargin==3)
    K=2;
end

%% Step 1. Create the expanded matrix with the new image
Aexp = [Uin*Din A]; % Column expansion
[Up,Dp] = qr(Aexp,0); % QR decomposition

%% Step 2. Compute SVD of Dp1
[Up1,Dp1,~] = svd(Dp,'econ');

%%  Step 3 Get the K greater values of Dp
D = Dp1(1:K,1:K);
U = Up*Up1(:,1:K);
