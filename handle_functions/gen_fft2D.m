function [F_tmp,varargout] = gen_fft2D(coef,tx,ty,F,varargin)


[F_tmp,fx] = gen_fft2(coef,tx,F(:,1,1,1,1,1),0); F0 = squeeze(F(1,:,1,1,1,1,1,1)).';
[F_tmp,fy] = gen_fft2(coef,ty,F0,0);

Nx = length(tx); Ny = length(ty);
num_arg = length(varargin);
if num_arg >= 1; fx = varargin{1}; end
if num_arg >= 2; fy = varargin{2}; end
if num_arg >= 3; Nx = varargin{3}; end
if num_arg >= 4; Ny = varargin{4}; end

vec_ind = 1:length(size(F)); 
vec_ind(1:2) = [2,1];

[F_tmp,fx] = gen_fft2(coef,tx,F,fx,Nx);
F_tmp = permute(F_tmp,vec_ind);
[F_tmp,fy] = gen_fft2(coef,ty,F_tmp,fy,Ny);

F_tmp = permute(F_tmp,vec_ind);

varargout{1} = fx(:); 
varargout{2} = fy(:).'; 

end


% Example 
% x = linspace(-7,7,101).';
% y = linspace(-7,6,201);
% F = exp(-x.^2 - (y-1).^2 ).*(x+2.*y);
% [G,fx,fy] = gen_fft2D(+1,x,y,F);
% [F2] = gen_fft2D(-1,fx,fy,G,x,y);
% fmx = [fx(1),fx(end)]; fmy = [fy(1),fy(end)];
% [GG,fx2,fy2] = gen_fft2D(+1,x,y,F,fmx,fmy,105,2001);
% [FF] = gen_fft2D(-1,fx2,fy2,GG,x,y);
% G_real = exp(-pi.*(2.*1i.* fy + pi .*(fx.^2 + fy.^2))).*pi.*(2 - 1i.*pi.*(fx + 2.* fy));
% 
% err_n(FF,F)
% err_n(G_real,G)








