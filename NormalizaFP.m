function [IN,M] = NormalizaFP(I,R)

% NORMALIZAFP performs the normalization of a fringe pattern using the
% technique described in Quiroga et al "Algorithm for fringe pattern
% normalization", Opt. Comm., 197, 2001
%
% @ 2021 Infor AOCG-UCM

% Initial operations. Set R value
if(nargin==1)
    %R = 2.5/max(size(I)); % Set R to a fixed value of 5% the Nyquist frequency
    R = 0.01;
end

% Create filters
tam = size(I);
u = linspace(-0.5,0.5,tam(2)); % Vector of horizontal frequencies
v = linspace(-0.5,0.5,tam(1)); % Vector of vertical frequencies
[U,V] = meshgrid(u,v); 
Masku = U>0;
Masku(1:round(tam(1)/2),round(tam(2)/2)+1) = 0;
Maskv = V>0;
Maskv(round(tam(1)/2)+1,1:round(tam(2)/2)) = 0;
MaskR = (U.^2+V.^2)>R.^2;
Masku = Masku.*MaskR;
Maskv = Maskv.*MaskR;

% Perform operations in frequency space
% Hamming window for suppressing ringing noise
HW = hamming(tam(1))*hamming(tam(2))';
F = fftshift(fft2(I)); % Fourier transform of intensity
FI1 = fftshift(F.*Masku.*HW);
FI2 = fftshift(F.*Maskv.*HW);
I1 = ifft2(FI1);
I2 = ifft2(FI2);

% Get the normalized fringe pattern and modulation
IN = (real(I1) + real(I2))./(abs(I1)+abs(I2));
M = abs(abs(I1) + 1i*abs(I2));