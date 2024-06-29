function y = heaviside2(x,varargin)

num_var = length(varargin); v0 = 0.5;
if num_var == 1; v0 = varargin{1}; end

y = x.*0;
y(x>0) = 1;

y( not(x == conj(x))) = NaN;

y(x ==0 ) = v0;


end

