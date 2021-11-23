function [varargout] = qcaDemod(varargin);

%QCADEMOD This function obtains the wrapped phase from a sequence of
%phase-shifted interferograms using the Quadrature Component Analysis
%algorithm (QCA). The mathematical details can be seen in [1]

% Usage: [pw,Mod,U1,U2,V] = pca(I,K,Mask)
% Inputs:
%   I  [NRows x NCols x num] Image 3D matrix where NRows and NCols are the 
%   number of rows and columns of the fringe patterns and num is the number 
%   of interferograms, I(:,:,1) is the first interfergram.   
%   K  1x1  Number of eigenvalues and eigenvectors to be returned.
%   Mask [NRows x NCols] is the processing mask
% Outputs:
%   pw  [NRows x NCols]  Wrapped modulating phase.
%   mod [NRows x NCols]  Modulation term.
%   Q1  [NRows x NCols]  First quadrature component.
%   Q2  [NRows x NCols]  Second quadrature component.

%   Javier Vargas, 
%   31/08/12 
%   Copyright 2012 
%   Biocomputing Unit, Centro Nacional de Biotecnología-CSIC 
%   $ Revision: 1.0.0.0 $
%   $ Date: 31/08/12 $

%NOTE: If you use or redistribute these functions, please refer the paper
%[1]
%
%REFERENCES
%
%[1] J. Vargas, C. O. S. Sorzano,
%  "Quadrature Component Analysis interferometry",
%   send for publication to Optics Letters

try
    funcName='qcaDemod';
    numvarargs = length(varargin);
    if numvarargs > 3
        error(['MyToolbox:' funcName ':TooManyInputs'], ...
            'requires at most 2 optional inputs: compound interferogram data and number of eigenvalues');
    end
    
    dim = size(varargin{1});
    %Compound image formed by the different interferograms columnwise
    %stacked
    X=reshape(varargin{1},dim(1)*dim(2),dim(3));
    %Number of eigenvalues and eigenvectors to be returned.
    K = varargin{2};
    Mask = varargin{3};
    
    for i=1:dim(3)
        temp = X(:,i);
        %Compound image after filtering the mask
        Xf(:,i)= temp(Mask(:));
    end
    
    [M N]=size(Xf);
    Xm=mean(Xf,2);
    Xd=Xf-repmat(Xm,1,N);
    C=Xd'*Xd;
    [a D at]=svd(C);    

    %Quadrature components
    U=Xd*a;
        
    U1 = X(:,1).*0;
    U2 = U1;    

    U1(Mask)=U(:,1);
    U2(Mask)=U(:,2);
    
    Q1 = reshape(U1,dim(1),dim(2));
    Q2 = reshape(U2,dim(1),dim(2));
    
    %modulating phase and modulation term
    pw = atan2(Q2,Q1);      
    Mod = sqrt(Q1.^2+Q2.^2);
    indPS = (((D(1,1)-D(3,3))/(D(2,2)-D(3,3))-1)/10)*100;
    indN = (D(3,3))/((D(1,1)+D(2,2)-2*D(3,3)))*100;
    ind =  sqrt( indPS.^2 +indN.^2);
    
    varargout{1} = pw;
    varargout{2} = Mod;
    varargout{3} = indPS;    
    varargout{4} = indN;
    varargout{5} = ind;    
            
catch ME
    throw(ME)
end
