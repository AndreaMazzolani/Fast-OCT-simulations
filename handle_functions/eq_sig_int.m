function [Ik,varargout] = eq_sig_int(Ik_raw,f_raw,varargin)
f_raw = f_raw(:); Nk = length(f_raw);


if nargin == 3;
    f_vec = varargin{1};
    Nkk = length(f_vec);
else
    Nkk = Nk;
    f_vec = linspace(f_raw(1),f_raw(Nk),Nkk);
end

if not(size(Ik_raw,1) == Nk)
    error('the first dimension of the field must equal the number of frequencies');
end
     
try 
    Ik = interp1(f_raw,Ik_raw,f_vec,'spline'); 
%     Ik = interp1(f_raw,Ik_raw,f_vec,'makima'); 
catch 
    Ik = interp1(f_raw,Ik_raw,f_vec,'pchip'); 

end
    varargout{1} = f_vec;

end