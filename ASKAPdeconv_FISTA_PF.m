function [Model Residual]=ASKAPdeconv_FISTA_PF(Dirtymap,PSF,center,lambda,niter,positiveflg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA -- a L1 norm based algorithom for solving the deconvolution problem in
%radio astronomy
% Modified on the 13th Sep 2010
% Details can be found in paper "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems"
%
% Input:
%
% 	Dirtymap the blurred image
%
% 	PSF ....    Point Spread Function of the Dirymap,here the psf is surposed to
%     	                be the same size of Dirtymap
%                            
% 	center      A vector indicating the peak coordinate of the PSF for
%		      example [129 129]
%                                         
% 	lambda    Regularization parameter, and also the noise level. However, if
%                        IUWT is adopted, the lambda should be around 10 times larger than the
%                         noise level, because the wavelet coefficients are scaled in the IUWT
% 
% 	niter       the number of iterations
%
% 	positiveflg  the flag of positive prior. If positiveflg=1, the source are
% 		          all positive in the model, by default positiveflg=0;
%
% Return:
%
% 	Model           the output cleaned image
%
% 	Residual       the residual image = Dirtymap-Model*PSF;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weight=sum(sum(PSF));
PSF=PSF/weight;

wh=size(Dirtymap) == size(PSF);
if wh(1) & wh(2)
else
  error(' The dirtymap and the dirty beam have to be the same size');
end

if nargin <=5
    positiveflg=0; 
end

%[m,n]=size(Dirtymap);

% computing the UV mask with the psf         
UV=fft2(circshift(PSF,1-center)); 
Fdirty=fft2(Dirtymap);

% initialization
X_temp=Dirtymap;
X=X_temp;
t_new=1;
for i=1:niter
    X_old=X_temp;
    t_old=t_new;
    
	% Gradient   
    D=UV.*fft2(X)-Fdirty;
    X=real(X-2/1*ifft2(conj(UV).*D));
    
    % Soft thresholding 
    Truncted=abs(X)-lambda;         
    index=find(Truncted<0);  % Find the non-important signals
    WShrink=sign(X).*((Truncted>0).*Truncted);
    WShrink(index)=0;        % Set the non-important signals to zero by force to make sure the backgroud of the model image is clean
    X_temp=WShrink; 

    %updating t and X
    t_new=(1+sqrt(1+4*t_old^2))/2;
    X=X_temp+(t_old-1)/t_new*(X_temp-X_old);

    Residual=Dirtymap-real(ifft2(UV.*fft2(X_temp)));
    %likelyhood=norm(Residual,'fro')^2;
    %total=likelyhood+lambda*sum(sum(abs((X_temp))));
    %fprintf('%3d    %15.5f  %15.5f  %15.5f    \n',i,norm(X_temp-X_old,'fro')/norm(X_old,'fro'),likelyhood,total);   
    
    if (positiveflg)
	   fprintf('Positive Flag Set');
       X_temp=X_temp.*(X_temp>0);
       X=X.*(X>0);
    end
	
end
Model=X_temp/weight ; 