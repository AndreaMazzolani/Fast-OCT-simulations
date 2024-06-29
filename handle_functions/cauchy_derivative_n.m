%% function calculating the n-th derivative using Cauchy formula
function [dyn] = cauchy_derivative_n(F,T,r,n);
T = T(:);
sF = size(F); Nt = length(T);
if (sF(1) ~= Nt)
    error('First dimension of the function f must be equal to the length of T'); 
end

Fs = exp(log(F) -2.*pi.*1i.*n.*T-n.*log(r));% F.*exp(-2.*pi.*1i.*n.*T)./r.^n; 
sn = sum(Fs(2:Nt-1,:,:,:,:),1)+1/2*(Fs(1,:,:,:,:)+Fs(Nt,:,:,:,:));
dyn = exp(log(sn)-log(Nt)+sum(log(1:n))); %sn/Nt*factorial(n);
end


