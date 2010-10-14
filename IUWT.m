function c = IUWT(x, nlevel)
%IUWT Isotropic Undecimated Wavelet Transform by Feng Li
%Reference paper: 
%"The Undecimated Wavelet Decomposition and its Reconstruction."
%Jean-Luc Starck, Jalal Fadili, and Fionn Murtagh
%  
 % Example:
%   Y = IUWT(X, 5);    % Decompose image X up to 5 level
%   R = IUWT(Y, -5);   % Reconstruct from Y
if nlevel>6   
   nlevel=6; 
end
    
 lp=[1/16,4/16,6/16,4/16,1/16	];  %Low band filter
%hp=[-1/16,-4/16,10/16,-4/16,-1/16]; %High band filter

bandnorm=[0.889434,  0.200105, 0.0857724, 0.0413447, 0.0202689, 0.00995628, 0.00513504];
% bandnorm=ones(1,10);
[m n]=size(x);


 
%----------------  remain unchanged when nlevel = 0  -------------------%
if nlevel == 0
    c = x;
%--------------------  decomposition,  if nlevel < 0  ------------------%
elseif nlevel > 0
     
    c = zeros(m,(nlevel+1)*n); %%simplified
    x = double(x);

    for i = 1:nlevel
        
        step=2^(i-1)-1;
        insert=zeros(1,step);
        lp=[1/16  insert  4/16 insert 6/16 insert 4/16 insert 1/16	];
        temp = symconv2(x, lp, 'col');    % low filtering

        ll = symconv2(temp, lp, 'row');   % low filtering
         
        highfrequen=x-ll;
        
        if i==nlevel
           ll=ll/bandnorm(i+1); 
        end
        
        c(1:m, (i-1)*n+1:(i+1)*n) = [highfrequen/bandnorm(i)  ll ]; %%simplified
       % replace x with ll for next level 
        x = ll;
        clear lp
        
    end
%--------------------  reconstruction,  if nlevel < 0  -----------------%
else
    nl = -nlevel;
    smln=n/(nl+1);%%simplified    
    for i = nl:-1:1
 
        ll = x(1:m, i*smln+1:(i+1)*smln);
        HVD=x(1:m, (i-1)*smln+1:i*smln);
        if i==nl
           ll=ll*bandnorm(i+1); 
        end
        x(1:m, (i-1)*smln+1:i*smln)= ll+HVD*bandnorm(i);
    end    
    % output
    c = x(:,1:smln);
end

%------------------------- internal function  --------------------------%
%       2-dimension convolution with edges symmetrically extended       %
%-----------------------------------------------------------------------%
function y = symconv2(x, h, direction)
 
l = length(h); s = size(x);
lext = (l-1)/2; % length of h is odd 
h = h(:)'; % make sure h is row vector 
y = x;
if strcmp(direction, 'row') % convolving along rows
    if ~isempty(x) && s(2) > 1 % unit length array skip convolution stage
        for i = 1: lext
            x = [x(:, 2*i), x, x(:, s(2)-1)]; % symmetric extension
        end
        x = conv2(x, h);
        y = x(:, l:s(2)+l-1); 
    end
elseif strcmp(direction, 'col') % convolving along columns
    if ~isempty(x) && s(1) > 1 % unit length array skip convolution stage
        for i = 1: lext 
            x = [x(2*i, :); x; x(s(1)-1, :)]; % symmetric extension
        end
        x = conv2(x, h');
        y = x(l:s(1)+l-1, :);
    end
end    
 