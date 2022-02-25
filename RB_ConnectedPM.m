function Wphi = RB_ConnectedPM(I)

% RB_CONNECEDPM Function for demodulating phase maps using the connected interferogram
% technique sescribed in Deng, et al, "Precise phase retrieval under harsh
% conditions..." Sci Rep, 2016
%
% Wphi= Rb_ConnectedPM(I) demoduates the phase by computing a set of new interferogram from
% the sums and differences of the original images stored in the 3D array I
%
% @ 2022, Infor, AOCG-UCM

% Initial operations
[R,C,N] = size(I);

% Computing sums and differences performing normalization
A = zeros(R,C,N*(N-1)/2);
S = zeros(R,C,N*(N-1)/2);
ind = 1;
for k=1:N
    for j=(k+1):N
        S(:,:,ind) = NormalizaFP(I(:,:,k)-I(:,:,j));
        A(:,:,ind) = NormalizaFP(I(:,:,k)+I(:,:,j));
        ind = ind+1;
    end
end

% Constructing Z
Z = cat(3,-S,A,S,-A);

% PCA algorithm
% Constructing X matrix
tam = size(Z);
X = zeros(tam(1)*tam(2),tam(3));
for k=1:tam(3)
    Z1 = Z(:,:,k);
    X(:,k) = Z1(:);
end
% Computing PCA
[~, score] = pca(X);
% Computing wrapped phase
Wphi = atan2(-score(:,2),score(:,1));
Wphi = reshape(Wphi,tam(1),tam(2));