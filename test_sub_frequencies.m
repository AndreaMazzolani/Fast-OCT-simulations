
clear all ;

desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus

if isfolder(desktop)
    path_str = desktop;
elseif isfolder(remote)
    path_str = remote
end

addpath(genpath([path_str,'handle_functions']));
addpath(genpath([path_str,'Illumination_functions']));


lam_min = 540.5e-9;  lam_max = 1007.8e-9;  % maximum wavelength (meters)
df_vec = 4600; % This is a choice that euristically works
df_long = 180; % This is a choice that euristically works
f_min = 1./lam_max; f_max = 1./lam_min;  % maximum frequency of the simulated spectrum
Nk = 7; ceil((f_max-f_min)./df_vec); Nlong = ceil((f_max-f_min)./df_long);
zmax1 = 20e-7;  Nk_sub = 4;
param = zeros(Nk,Nlong);
param_sub = zeros(Nk_sub,Nlong);
f_tot_sub = zeros(Nk_sub,Nlong);
ind_tot_sub = zeros(Nk_sub,Nlong);
err_tot = zeros(Nlong,1); err_tot_sub = zeros(Nlong,1);

df = (f_max-f_min)/(Nk-1); espa = 1;
f_vec = f_min + (linspace(0,(f_max-f_min).^(1/espa),Nk)).^espa;

f_long = linspace(f_min ,f_max,Nlong);
f0 = (f_min + f_max)/2;
Temps = 0;   zmax2 = 44e-6; a = (2*pi*zmax1); error = 1e-9 ; no_itr = 20584 ;

%  NLONG = flip(1:Nlong);
itt = 0; ind_sub = 1:Nk_sub; testf1 = 0;

if0 = ceil(Nlong*sqrt(2)/2);
v = f_long(if0);

if_short = max(find(f_vec <= v));
if isempty(if_short); if_short = 1; end;
xm = zeros(2.*Nk,1);

if if_short == Nk
    xm(if_short) = 1;
else
    vm = f_vec(if_short); vM = f_vec(if_short+1);
    xm(if_short) = (vM-v)./(vM-vm);  xm(if_short+1) = (v-vm)./(vM-vm);
end

if_slit = 0;
if if_short + Nk_sub-1 > Nk
    if_slit = if_short + Nk_sub-1 - Nk;
end

%%%%%%%%%%%%%%%%     ERROR FUNCTION TO MINIMIZE:   %%%%%%%%%%%%%%%%%%%%
ink2 = 0; %if (if_short > Nk_sub) & (if_short < Nk - Nk_sub); ink2 = floor((Nk_sub-1)/2); end
if (if_short > Nk_sub) & (if_short < Nk - Nk_sub); ink2 = floor((Nk_sub-1)/2); end
Ntt = [Nk,Nk_sub];
ind_sub_old = ind_sub;
ind_sub = [(if_short -ink2): (if_short + Nk_sub-1-ink2)]- if_slit;

if not(err_n(ind_sub,ind_sub_old) == 0)
    itt = itt+1;
    indt(itt) = if0;
end
f_sub = f_vec(ind_sub);
ind_tot_sub(:,if0) = ind_sub;

fn = @(x) loss_freq(x,a,v,f_vec); Dfn = @(x) der_loss_freq(x,a,v,f_vec);
fn_sub = @(x) loss_freq(x,a,v,f_sub); Dfn_sub = @(x) der_loss_freq(x,a,v,f_sub);

%%%%%%%%%%%%%%%%%%%%%%     OPTIMISATION     %%%%%%%%%%%%%%%%%%%%%%


xm_sub = zeros(2*Nk_sub,1); xm_sub(1:Nk_sub) = xm(ind_sub); xm_sub(Nk_sub+1:2*Nk_sub) = xm(Nk+ind_sub);

[v2 , no_itr22, norm1,lambda] = grad_desc(xm,fn,Dfn,no_itr,error*1e1);
[v2_sub , no_itr_sub, norm1_sub,lambda_sub] = grad_desc(xm_sub, fn_sub, Dfn_sub, no_itr*100,error);

 
for ik = 1:Nk
    param(ik,if0) = v2(ik).*exp(1i.*v2(Nk+ik));
    %             lamb_opt(ik,if0) = x_opt(ik).*exp(1i.*x_opt(Nk+ik));
end

for ik = 1:Nk_sub
    param_sub(ik,if0) = v2_sub(ik).*exp(1i.*v2_sub(Nk_sub+ik));
end
f_tot_sub(:,if0) = f_sub;

param2 = param(:,if0); param2_sub = param_sub(:,if0); %fL= FF.f_long; v = fL(31);
%         lamb2_opt = lamb_opt(:,if0); %fL= FF.f_long; v = fL(31);

N0 = 1501; A = 0;rand(N0,1); zm = zmax1; B = zm.* sort(rand(N0,1))*1;
Sk = exp(-(v./f0).^2.*A.^2).*exp(1i.*(B).*2.*pi.*v);
S_est = sum(param2.'.*exp(-(f_vec./f0).^2.*A.^2).*exp(1i.*(B).*2.*pi.*f_vec),2);
S_est_sub = sum(param2_sub.'.*exp(-(f_sub./f0).^2.*A.^2).*exp(1i.*(B).*2.*pi.*f_sub),2);

comp(S_est_sub,Sk,'',B*1e6)
err_n(S_est_sub,Sk)

% figure(1); clf;1./f_sub(1:4)*1e9

err_tot(if0) = err_n(Sk,S_est);
err_tot_sub(if0) = err_n(Sk,S_est_sub);
E_k = exp(1i.*(B).*2.*pi.*v);
E_k_est = sum(param2_sub.'.*exp(1i.*(B).*2.*pi.*f_sub),2);
E_sub = zeros(N0,Nk_sub); st_s = 'E_k_est(:)';
for ik_s = 1:Nk_sub
    iks = num2str(ik_s);
    eval(['E_sub(:,',iks,') = exp(1i.*(B).*2.*pi.*f_sub(',iks,'));']);
    st_s = [st_s,', E_sub(:,',iks,')'];
end
%        st_s = [st_s,'1'];
% eval(['comp(',st_s,')'])
% 
% comp(E_sub(:,1),E_k)



figure(1); clf; NF = 12; NL = 1.4;
subplot(3,2,1); plot(B*1e6,real(E_sub(:,1)),'LineWidth',NL);
subplot(3,2,3); plot(B*1e6,real(E_sub(:,2)),'LineWidth',NL);
subplot(3,2,5); plot(B*1e6,real(E_sub(:,3)),'LineWidth',NL);
subplot(3,2,6); plot(B*1e6,real(E_sub(:,4)),'LineWidth',NL);
subplot(3,2,2); plot(B*1e6,real(E_k_est),'k','LineWidth',NL); hold on; plot(B*1e6,real(E_k),'--r','LineWidth',NL*1.8); legend('Rigorous','Estimated');




