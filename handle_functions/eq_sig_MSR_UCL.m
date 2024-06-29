function [Ik,varargout] = eq_sig_MSR_UCL(Ik_UCL,varargin)

c = 2.997924580105029e+08; 
path_MSR = 'OCT_models/Simulated_OCT_DW_UL_Bscan/';
addpath(genpath('handle_functions'));

% DATAF = load([path_MSR,'OCT_raw_eq_spec_Nk_2048_Nksub_5_Nlong_7580_zm_2.5e-05_lam_min_1170_lam_max_1408.mat']);  
% DATAF = load([path_MSR,'OCT_raw_eq_spec_Nk_2048_Nksub_11_Nlong_9601_zm_0.0005_lam_min_1170_lam_max_1408.mat']);  
%  DATAF = load([path_MSR,'OCT_raw_eq_spec_Nk_2048_Nksub_7_Nlong_9601_zm_0.0001_lam_min_1170_lam_max_1408.mat']);  
% DATAF = load([path_MSR,'OCT_raw_eq_spec_Nk_2048_Nksub_31_Nlong_9601_zm_0.0015_lam_min_1170_lam_max_1408.mat']);  
DATAF = load([path_MSR,'OCT_raw_eq_spec_Nk_2048_Nksub_151_Nlong_6858_zm_0.0055354_lam_min_1170_lam_max_1408.mat']);

if isfield(DATAF,"DATAF"); DATAF = DATAF.DATAF; end; 
% unpackStruct(DATAF) 

f_raw = DATAF.f_vec; Nk_raw = length(f_raw);
% f_raw = f_raw_tmp(ind_UCL);

if nargin == 2;
    f_vec = varargin{1};
    Nk = length(f_vec);
else
    Nk  = Nk_raw;
    f_vec = linspace(f_raw(1),f_raw(Nk_raw),Nkk);
end

if not(size(Ik_UCL,1) == Nk_raw)
    error('the first dimension of the field must equal the number of frequencies');
end

f_min = min(f_vec); f_max = max(f_vec);
f_long_tmp = DATAF.f_long; ind_min_long = max(find(f_long_tmp <= f_min)); ind_max_long = min(find(f_long_tmp >= f_max)); IND_long = ind_min_long:ind_max_long;
f_short_tmp = DATAF.f_vec; Nk_short_tmp = length(f_short_tmp); 
ind_min_short = max(1,max(find(f_short_tmp <= f_min))-ceil(Nk_short_tmp./2)-1); ind_max_short = min(Nk_short_tmp,min(find(f_short_tmp >= f_max))+ceil(Nk_short_tmp./2)+1); IND_short = ind_min_short:ind_max_short;

lamb_low = 1e9./f_long_tmp(end); lamb_high = 1e9./f_long_tmp(1);
if (isempty(ind_min_long)||isempty(ind_max_long)||isempty(ind_min_short)||isempty(ind_max_short))
    error(['The spectrum is too large! The available interval of wavelengths is [,',num2str(lamb_low,5),',',num2str(lamb_high,5),'] nanometers. Run another suitable MSR!']);
end

f_long = f_long_tmp(IND_long); Nk_long = length(f_long);
f_short = f_short_tmp(IND_short); Nk_short = length(f_short);   
f_long_t = f_long.*c; f_short_t = f_short.*c;
param = zeros(Nk_short,Nk);
param_tot_tmp = DATAF.param; param_tot = param_tot_tmp(IND_short,IND_long);
param_tot_sub_tmp = DATAF.param_sub; param_tot_sub = param_tot_sub_tmp(:,IND_long);
param_tot_subt = DATAF.param_subt;
Nk_sub = size(param_tot_sub,1); param_sub = zeros(Nk_sub,Nk);
f_tot_sub_tmp = DATAF.f_tot_sub; f_tot_sub = f_tot_sub_tmp(:,IND_long);
ind_tot_sub_tmp = DATAF.ind_tot_sub; ind_tot_sub = ind_tot_sub_tmp(:,IND_long);
ind_tot_subt = DATAF.ind_tot_subt;
% ind_tot_subt_tmp = DATAF.ind_tot_subt; ind_tot_subt = ind_tot_subt_tmp(:,IND_long);

% amp_par = (f_vec(:).')./f_short(:);

indt = DATAF.indt; %indt(1) = 1; indt(end) =  size(param_tot_sub,2);;
Nindt = length(indt);


for ik_short = 1:Nk_short
    param(ik_short,:) = interp1(f_long,param_tot(ik_short,:),f_vec);
end

if sum(isnan(param(:))) > 0
    error(['The spectrum is too large! The interval of wavelength available is [,',num2str(lamb_low,5),',',num2str(lamb_high,5),'] nanometers']);
end

% Piecewise f_sub interpolation %%%  ind_tot_sub = zeros(Nk_sub,Nlong);
itt = 1; f_sub_vec = zeros(Nk_sub,Nk);  ind_sub_vec = zeros(Nk_sub,Nk);
for ikL = 1:(Nk_long-1)
    fL1 = f_long(ikL); fL2 = f_long(ikL+1);
    indL_k =  find((f_vec >= fL1)&(f_vec <= fL2));
    ind_sub1 = ind_tot_sub(:,ikL); ind_sub2 = ind_tot_sub(:,ikL+1);
    if not(isempty(indL_k))
        if ind_sub1(1) ==  ind_sub2(1) % same base f_sub
            for ik_sub = 1:Nk_sub
                param_sub(ik_sub,indL_k) = interp1(f_long(ikL:(ikL+1)),param_tot_sub(ik_sub,ikL:(ikL+1)),f_vec(indL_k)); % interpolation over each fixed "short frequency"
            end
            %             f_sub_vec(:,indL_k) = f_tot_sub(:,ikL).*ones(Nk_sub,length(indL_k));
            %             ind_sub_vec(:,indL_k) = ind_sub1(:).*ones(Nk_sub,length(indL_k));

        else % DIFFERENT base f_sub, I have to use the second calculation for the right-hand border

            itt = min(find(indt > IND_long(ikL))); ind_sub2 = ind_tot_subt(:,itt);
                        if not(ind_sub1(1) ==  ind_sub2(1));
                            error(['Array of indices is made wrongly!']);
                        end
            for ik_sub = 1:Nk_sub
                param_sub(ik_sub,indL_k) = interp1(f_long(ikL:(ikL+1)),[param_tot_sub(ik_sub,ikL),param_tot_subt(ik_sub,itt)],f_vec(indL_k)); % interpolation over each fixed "short frequency"
            end
            %             itt = itt +1;
        end
        f_sub_vec(:,indL_k) = f_tot_sub(:,ikL).*ones(Nk_sub,length(indL_k));
        ind_sub_vec(:,indL_k) = ind_sub1(:).*ones(Nk_sub,length(indL_k));
    end
end




Ik_sub = Ik_UCL(ind_sub_vec);


 Ik = sum((param_sub).*Ik_sub);


    varargout{1} = f_vec;

end