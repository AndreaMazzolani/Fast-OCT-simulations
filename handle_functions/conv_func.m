
function [fg_c,y] = conv_func(x,f,g,varargin);

num_arg = length(varargin);
Nx = length(x); [F_tmp,vecy_tmp] = gen_fft2(1,x(:,1,1,1,1,1,1),f(:,1,1,1,1,1,1),0);
vecy = [vecy_tmp(1),vecy_tmp(Nx)]; Ny = Nx;
if num_arg >= 1; vecy = varargin{1}; end
if num_arg >= 2; Ny = varargin{2}; end

% x = x(:); f = f(:); g = g(:);
[F,y] = gen_fft2(1,x,f,vecy,Ny);
[G,y] = gen_fft2(1,x,g,vecy,Ny);
fg_c = gen_fft2(-1,y,F.*G,x);


end



%%% TEST 
 
% x = linspace(-5,5,201).';
% f = exp(-2.*abs(x)); g = exp(-(x-1).^2);
% fg_c_anal = sqrt(pi)./2.*exp(-2.*x-1).*(exp(4).*erfc(2-x) + exp(4.*x).*erfc(x));
% comp(f,g,fg_c_anal,'',x)
% vecy = [-5,5]; Ny = 201;
% 
% [fg_c,y] = conv_func(x,f,g,vecy,Ny);
% 
% comp(fg_c,fg_c_anal,'',x)
% 
% err_n(fg_c,fg_c_anal)
