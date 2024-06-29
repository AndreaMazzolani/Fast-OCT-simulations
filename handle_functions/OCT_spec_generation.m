function [Sk_raw] = OCT_spec_generation(f_vec , varargin)

Ng = 5;
if nargin == 2
    Ng = varargin{1}; if Ng <= 0; error('The second input must be a poistive integer number'); end
end

f_vec = f_vec(:);
f_max =  max(f_vec); f_min =  min(f_vec);  Nk = length(f_vec); f0 = (f_min + f_max)./2; cut_off = 1e-3;  Wk = 2*pi*(f_max-f0)/sqrt(log(1/cut_off));

DFF0 =  (f_max-f_min)*0.1; DF0 = f_max-f_min-2*DFF0; Noise_Spf = 0.1.*(1+ randn(Nk,Ng))/Ng;  rg = 1+ 0.2*(2*rand(1,Ng)-1); expg = rand(1,Ng)*1.4+1; Wrg = rand(1,Ng);
Wkg = Wk*(1+Wrg); f0g_norm = linspace(0,1,Ng)  + 0.1*(2*sort(rand(1,Ng))-1); f0g = f_min + DFF0 + DF0.*f0g_norm;
SK_raw = rg.*exp(-((2*pi./Wkg).*abs(f_vec-f0g)).^expg).*(1+Noise_Spf);
Sk_raw = sum(SK_raw,2);

% clf; plot(f_vec,Sk_raw)
end