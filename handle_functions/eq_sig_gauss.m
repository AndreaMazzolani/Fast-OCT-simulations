function [Ik,varargout] = eq_sig_gauss(Ik_raw,f_raw,varargin)

DATAF = load('zeus_F_RECR_opt_Nk_40_Nlong_1001_zm_4e-05.mat');
Nk = length(f_raw);  f_min = f_raw(1);  f_max = f_raw(Nk);
f_short_tmp = DATAF.f_vec; Nk_short_tmp = length(f_short_tmp);
ind_min_short = max(1,max(find(f_short_tmp <= f_min))-ceil(Nk_short_tmp./2)-1);
ind_max_short = min(Nk_short_tmp,min(find(f_short_tmp >= f_max))+ceil(Nk_short_tmp./2)+1);
IND_short = ind_min_short:ind_max_short;

f_long_tmp = DATAF.f_long; ind_min_long = max(find(f_long_tmp <= f_min)); ind_max_long = min(find(f_long_tmp >= f_max)); IND_long = ind_min_long:ind_max_long;
f_long = f_long_tmp(IND_long); Nk_long = length(f_long); f_short = f_short_tmp(IND_short); Nk_short = length(f_short);
param_tot_tmp = DATAF.param; param_tot = param_tot_tmp(IND_short,IND_long);

param = zeros(Nk_short,Nk);

for ik_short = 1:Nk_short
    param(ik_short,:) = interp1(f_long,param_tot(ik_short,:),f_raw);
end

f_raw = f_raw(:); Nk = length(f_raw);
Nkk = Nk;

f_vec = linspace(f_raw(1),f_raw(Nk),Nkk);

if not(size(Ik_raw,1) == Nk)
    error('the first dimension of the field must equal the number of frequencies');
end

% Ik = interp1(f_raw,Ik_raw,f_vec,'cubic');
sk = size(Ik_raw); Lsk = length(sk);
Ik = zeros(size(Ik_raw));
if Lsk <= 2; Ik =  (param.')*Ik_raw;
else
    if Lsk == 3
    for isk1 = 1:sk(2)
        for isk3 = 1:sk(3)
            Ik(:,:,isk2) = (param.')*Ik_raw(:,:,isk3);
        end
    end
    elseif Lsk == 4
        for isk3 = 1:sk(3)
            for isk4 = 1:sk(4)
                 Ik(:,:,isk3,isk4) = (param.')*squeeze(Ik_raw(:,:,isk3,isk4));
            end
        end
    end
end



varargout{1} = f_vec;


end