function [Sk] = spec_estim(f_vec,Ik)


% This estimation works only if the scatterers are not placed close to the mirror position.

df = diff(f_vec); df2 = diff(df); size_f = size(f_vec); sk = size(Ik);

if prod(size_f) == max(size_f)
    Nk = length(f_vec); f_vec = f_vec(:);
else
    error('The first input variable must be a one-dimensional array of frequencies')
end

if not(sk(1) == Nk); 
    error('The first dimension of the Dataset must be related to the frequencies');
end

Df = f_vec(Nk)-f_vec(1); 
Az0 = abs(ff_tr(Ik,f_vec,0,-1));  min_Az0 = mean(Az0);  global Azz; Azz = Az0;
Ikn = Ik./Az0.*min_Az0;  Sk_tmp = mean(Ikn,2);  Sk = Sk_tmp;
Nslit = 1;ceil(Nk/1500)*2+1; 
[Spec_test ,is] = slit_mean(Sk_tmp,Nslit); Sk = Sk_tmp; Sk(is) = Spec_test;
end



