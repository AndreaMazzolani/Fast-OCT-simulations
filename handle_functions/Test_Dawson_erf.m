clear all; clc; 


A = 13*(rand(10) + 1i.*rand(10));

Ea = double(erf(sym(A)));
Ea1 = erf_app(A);

Ea_D = -2.*1i./sqrt(pi).*dawson(1i.*A).*exp(-A.^2); % Erf as a function of dawson


err_n(Ea_D,Ea)
err_n(Ea1,Ea)
