function [a,b,c,varargout] = quad_interp(x,y,varargin)
% solution of min{||y- a+bx+cx^2||} , output = a,b,c. The regression is performed along the first dimension
sx = size(x); sy = size(y); W = ones(size(y));

Nx = sy(1); if (length(sx) == 2)  & (sx(1) == 1); x = x(:);  end
x = x.*ones(size(y));
% if not(isequal(sx,sy))
%     if (sy(1) == Nx); x = x.*ones(size(y)); else; error('the size of "x" and "y" arrays must match!'); end;
%     else; error('the size of "x" and "y" arrays must match!');
%     end
% end

num_var = length(varargin); 
if num_var == 1; W = varargin{1}; sW = size(W); if not(isequal(sW,sy)); error('the size of "x" and "W" arrays must match!'); end; end
if num_var > 1;  error('Too many input data'); end



if min(W(:))<0; error('Weights must be positive'); end;
W = W./sum(W);
Ex = sum(W.*x); Ex2 = sum(W.*x.^2); Ex3 = sum(W.*x.^3); Ex4 = sum(W.*x.^4); Ey = sum(W.*y); Eyx = sum(W.*y.*x); Eyx2 = sum(W.*y.*(x.^2));
S11 = (Ex2- Ex.^2);
S12 = (Ex3-Ex.*Ex2);
S22 = (Ex4-Ex2.^2);
Sy1 = (Eyx-Ey.*Ex);
Sy2 = (Eyx2-Ey.*Ex2);

b = (Sy1.*S22-Sy2.*S12)./(S22.*S11-(S12.^2));
c = (Sy2.*S11-Sy1.*S12)./(S22.*S11-(S12.^2));
a = Ey-b.*Ex-c.*Ex2;

% if sum(isnan(a(:)))+sum(isnan(b(:)))+sum(isnan(c(:))) > 0; error('There are NaN values, the regression failed!'); end

if nargout >= 4; varargout{1} = a + b.*x + c.*(x.^2); end
if nargout >= 5; error('The maximum variables in output must be 2'); end

end



%%%%%%% TEST
% % Example no weights
% x = linspace(-4,4,1001).';
% A0 = rand(1,12,3,4); B0 = rand(1,12,3,4); C0 = rand(1,12,3,4);
% y = A0.*ones(size(x)) + B0.*x + C0.*x.^2;
% [a,b,c] = quad_interp(x,y);
% 
% [err_n(a,A0)
% err_n(b,B0)
% err_n(c,C0)]
% 
% %% Example with W = 1; 
% X = sort(x(1) + rand(size(y)).*10.*(x(end)-x(1)));
% Y =  A0.*ones(size(X)) + B0.*X + C0.*X.^2;
% [A,B,C,Y_est] = quad_interp(X,Y);
% 
% [err_n(A,A0) err_n(B,B0) err_n(C,C0), err_n(Y,Y_est)]
% 
% %% Example with W general 
% W = (Y).^2;
% Y =  A0.*ones(size(X)) + B0.*X + C0.*X.^2;
% [A,B,C,Y_est] = quad_interp(X,Y,W);
% 
% [err_n(A,A0) err_n(B,B0) err_n(C,C0), err_n(Y,Y_est)]
% 
% 
% 
