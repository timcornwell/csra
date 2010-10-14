function [Model Residual]=FISTA_PF_IUWT(Dirtymap,PSF, center,lambda,niter,  positiveflg, waveletflg,level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA -- a L1 norm based algorithom for solving the deconvolution problem in
%radio astronomy
% Modified on the 13th Sep 2010
% Details can be found in paper "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems"

% Dirtymap the blurred image
%
% PSF .... Point Spread Function of the Dirymap,here the psf is surposed to
%           be the same size of Dirtymap
%                            
% center      A vector indicating the peak coordinate of the PSF for
% example [129 129]
%                                         
% lambda  Regularization parameter, and also the noise level. 

% niter    the number of iterations
%
% positiveflg  the flag of positive prior. If positiveflg=1, the source are
% all positive in the model, by default positiveflg=0;

% waveletflg   the flag of Istropic Undecimated Wavelet Transform (IUWT). If
% waveflg=1, IUWT is adopted, otherwise the partial Fourier

% level   the level of the wavelet transform, for this case, the level is
% no larger than 6. By default, it is set to  6
 

% Model the output cleaned image
% Residual the residual image = Dirtymap-Model*PSF;
% 

   
weight=sum(sum(PSF));
PSF=PSF/weight;
wh=size(Dirtymap) == size(PSF);
if wh(1) & wh(2)
else
  error(' The dirtymap and the dirty beam have to be the same size');
end

if nargin <=6
    waveletflg=0;
    level=6; 
  
end

if nargin <=7
    level=6; 
 
end
 

if waveletflg
W = @(x) IUWT(x,level);% W means IUWT
WT = @(x) IUWT(x,-level);% WT  means inverse IUWT
else
W = @(x) x;
WT = @(x) x; 
end




[m,n]=size(Dirtymap);
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
    X=X-real(ifft2(conj(UV).*D));
    WX=W(X);
    
    
    % Soft thresholding 
    % For IUWT, we need to protect the Low frequency components when doing
    % the soft thresholding
     
    if waveletflg   %IUWT, assuming the source are sparse in the IUWT domain
         Highfre=WX(1:m, 1:(level)*n); % High frequency components
         Truncted=abs(Highfre)-lambda;
         WShrink=sign(Highfre).*((Truncted>0).*Truncted);
         index=find(Truncted<0);
         WShrink(index)=0;
         WShrink=[WShrink  WX(1:m, 1+level*n:(level+1)*n) ];
   
    else           % Partial Fourier, assuming the source are sparse in the image domain
         Truncted=abs(WX)-lambda;         
         index=find(Truncted<0);  % Find the non-important signals
         WShrink=sign(WX).*((Truncted>0).*Truncted);
         WShrink(index)=0;        % Set the non-important signals to zero by force to make sure the backgroud of the model image is clean
         
    end
    X_temp=WT(WShrink); 
    %updating t and X
    t_new=(1+sqrt(1+4*t_old^2))/2;
    X=X_temp+(t_old-1)/t_new*(X_temp-X_old);
    Residual=Dirtymap-real(ifft2(UV.*fft2(X_temp)));
    likelyhood=norm(Residual,'fro')^2;
    total=likelyhood+lambda*sum(sum(abs(W(X_temp))));
    fprintf('%3d    %15.5f  %15.5f  %15.5f  %15.5f    \n',i,norm(X_temp-X_old,'fro')/norm(X_old,'fro'),likelyhood,total, max(max(abs(Residual))));   
    
    if (positiveflg)
    X_temp=X_temp.*(X_temp>0);
    X=X.*(X>0);
    end
 
    
end
    X_temp=X_temp.*(X_temp>lambda);
   
Model=X_temp/weight; 

