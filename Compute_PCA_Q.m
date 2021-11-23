function gamma = Compute_PCA_Q(D)

% COMPUTE_PCA_Q Computes the quality parameter from the eingenvalues
% assocciated with the principal components
%
% gamma = Compute_PCA_Q(D) computes the quality parameter, gamma, from the
% eigenvalues matrix D associated with a PCA analysis according to the
% algorithm (equation 28) given in Vargas et al, "Error analysis of the
% principal component analysis demodulation algorithm", Applied Physics B,
% 115, 355-364, (2014)
%
% @ 2021, Jose A. Gomez-Pedrero, AOCG-UCM
%
% DISCLAIMER: Disclaimer: This software is supplied without any guarantee 
% or any support whatsoever.

% Checks input matrix
if(not(numel(D)==9))
    error('Input should be a 3x3 square matrix')
end
X = (D(1,1)/D(2,2)-1)/10;
Y = D(3,3)/(D(1,1)+D(2,2)-2*D(3,3));
gamma = abs(X+1i*Y);

end