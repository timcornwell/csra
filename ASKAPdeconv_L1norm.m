function [Model]=ASKAPdeconv_L1norm(Dirtymap,PSF,center,lambda,niter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA -- a L1 norm based algorithom for solving the deconvolution problem in
%radio astronomy
%
% Details can be found in paper "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems"
%
% Input:
%
% 	Dirtymap the blurred image
%
% 	PSF .... Point Spread Function of the Dirymap,here the psf is surposed to
%           	  be the same size of Dirtymap
%                            
% 	center      A vector indicating the peak coordinate of the PSF for
% 		      example [129 129]
%                                         
% 	lambda  Regularization parameter, the larger the less iterations
% 	              required, however, the smaller the finer the reconstruction. 
% 
% 	niter    the number of iterations
%
% Return:
%
%	 Model the output cleaned image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSF=PSF/sum(sum(PSF)); % 
wh=size(Dirtymap) == size(PSF);
if wh(1) & wh(2)
else
  error('The dirtymap and the dirty beam have to be the same size');
end

if nargin < 5
    niter=100;
end

%[m,n]=size(Dirtymap);

% computng the UV mask with the psf         
UV=fft2(circshift(PSF,1-center));
Fdirty=fft2(Dirtymap);

%Calculate the Lipschitz constant as introduce in the paper
L=2*max(max(abs(UV).^2));
fprintf('Lipschitz constant : %f15.5',L);
pause;

% initialization
X_temp=Dirtymap;
X=X_temp;
t_new=1;
for i=1:niter
    X_old=X_temp;
    t_old=t_new;
	
    % Gradient
    D=UV.*fft2(X)-Fdirty;
    X=real(X-2/L*ifft2(conj(UV).*D));
	
    % Soft thresholding 
    D=abs(X)-lambda/(L);
    X_temp=sign(X).*((D>0).*D);
	
    %updating t and X
    t_new=(1+sqrt(1+4*t_old^2))/2;
    X=X_temp+(t_old-1)/t_new*(X_temp-X_old);
    
end
Model=X_temp;
