function [Ex,varargout] = Ex_Debye_freq_rad_Taylor(v,t0,f_vec_t,refind,w0,P,K0,th,t_w)
 
global Nderz_glob Nderx_glob dz_glob dx_glob;
c = 2.997924580105029e+08;  Nv = size(v,1);
DATAF = load('zeus_F_RECR_opt_Nk_35_Nlong_10001_zm_1.5e-05.mat');
% DATAF = load('zeus_F_RECR_opt_Nk_40_Nlong_1001_zm_4e-05.mat');
%DATAF = load('zeus_F_RECR_opt_Nk_48_Nlong_16667_zm_1.5e-05_lam_min_1000_lam_max_1500.mat'); 
f_short = DATAF.f_vec; Nk_short = DATAF.Nk; f_short_t = f_short.*c;
f_long = DATAF.f_long; N_long = length(f_long);
f_vec = f_vec_t./c; Nk = length(f_vec);
param_tot = DATAF.param;  param = zeros(Nk_short,Nk);
amp_par = (f_vec(:).')./f_short(:);

for ik_short = 1:Nk_short
    param(ik_short,:) = interp1(f_long,param_tot(ik_short,:),f_vec);
end
iv = (nargout == 2);


lamb_low = 1e9./f_long(N_long);
lamb_high = 1e9./f_long(1);

if sum(isnan(param(:))) > 0
    error(['The spectrum is too large! The interval of wavelength available is [,',num2str(lamb_low,5),',',num2str(lamb_high,5),'] nanometers']);
end

if Nk <= Nk_short
    MSR = 0; % No - Multi Spectral Regression
    f_short = f_vec; f_short_t = f_vec_t; Nk_short = Nk;
else 
    MSR = 1; % Multi Spectral Regression
end

dr_taylor = 10e-6; dr = dr_taylor/sqrt(2);
dz = 20e-6; 
Nderz = 3; Nderx = 13;

epsz = 1e-9; r_source_max = max(sqrt(v(:,1).^2 + v(:,2).^2)); %;%min(max(sqrt(v(:,1).^2 + v(:,2).^2)),100e-6); %100e-6;  max(sqrt(v(:,1).^2 + v(:,2).^2));%
if not(isempty(Nderx_glob)); Nderx = Nderx_glob; end;
if not(isempty(Nderz_glob)); Nderz = Nderz_glob; end;
if not(isempty(dz_glob)); dz = dz_glob; end;
if not(isempty(dx_glob)); dr = dx_glob; end;

r_source = [dr:2*dr:r_source_max+dr];
Z_hat = 0; minZ = min(v(:,3)); maxZ = max(v(:,3));
z_source = [minZ - dz: 2*dz-epsz: maxZ + dz]; [R_S,Z_S] = meshgrid(r_source,z_source); 
Nv_S = numel(R_S);
v_source = zeros(Nv_S,3); v_source(:,1) = R_S(:); v_source(:,2) = 0.*R_S(:); v_source(:,3) = Z_S(:);

[EJ0_Txz,EJ2_Txz] = Ex_Debye_freq_rad_der_x_z(v_source,t0,f_short_t,refind,w0,P,K0,th,t_w,Z_hat,Nderz,Nderx);

% Taylor expansion calculation
x = v(:,1); y = v(:,2); zP = v(:,3);
rP = sqrt(x.^2 + y.^2); phiP = reshape(atan2(y,x),[1,Nv]); 
Ex_short = zeros(Nk_short,Nv); Ey_short = zeros(Nk_short,Nv);

for ivs = 1:Nv_S
%     disp([num2str(ivs),'/',num2str(Nv_S)])
    vsLx = v_source(ivs,1); vsLy = v_source(ivs,2) ; vsLz = v_source(ivs,3) ;  vsLrho = sqrt(vsLx.^2 + vsLy.^2);
    DX1 = rP - vsLrho;
    ind_scanL = find( (zP <  vsLz + dz ) & (zP  >=  vsLz - dz ) & (rP < vsLrho + dr)& (rP >= vsLrho - dr));

    if not(isempty(ind_scanL))
        NL = length(ind_scanL); DX = DX1(ind_scanL); DZ = vsLz - zP(ind_scanL);
        TDX3 = transpose((DX.^(0:Nderx)./factorial(0:Nderx)));
        TDZ3 = transpose(((Z_hat-DZ).^(0:Nderz)./factorial(0:Nderz)));
        Phase_freq_term = exp(1i.*refind.*2.*pi.*f_short(:).*( Z_hat.*ones(1,NL)-DZ.' ));
        Ex_vec_tmpj0 = zeros(Nk_short,NL);  Ex_vec_tmpj2 = Ex_vec_tmpj0;

        for iza2 = 1:Nderz+1
            Ex_vec_tmpj0 = Ex_vec_tmpj0 + (reshape(EJ0_Txz(:,ivs,iza2,:),[Nk_short,Nderx+1])*TDX3).*(TDZ3(iza2,:));
            Ex_vec_tmpj2 = Ex_vec_tmpj2 + (reshape(EJ2_Txz(:,ivs,iza2,:),[Nk_short,Nderx+1])*TDX3).*(TDZ3(iza2,:));
        end

        Ex_short(:,ind_scanL) = squeeze(Ex_vec_tmpj0 +  Ex_vec_tmpj2.*cos(2.*phiP(1,ind_scanL))).*Phase_freq_term;
        if iv 
            Ey_short(:,ind_scanL) = squeeze(Ex_vec_tmpj2.*sin(2.*phiP(1,ind_scanL))).*Phase_freq_term;
        end
    end
end

p2 = param.*amp_par;
zvec = reshape(zP,[1,Nv]);
Phase_short = exp(-1i.*2.*pi.*refind.*zvec.*f_short(:));
Phase_long = exp(1i.*2.*pi.*refind.*(zP.').*(f_vec(:)));

% Ex_vec_tmp_short3 = Ex_vec_tmp_short.*Phase_short;
% Ex_vec_tmp = squeeze(sum(Ex_vec_tmp_short3.*p2,1)).*Phase_long; toc;

if MSR
    Ex_short2 = reshape(Ex_short,[Nk_short,Nv]); Ex_short3 = Ex_short2.*Phase_short;
    Ex = zeros(Nk,Nv);
    for ik = 1: Nk 
        Ex(ik,:) = ((p2(:,ik)).')*(Ex_short3.*Phase_long(ik,:));
    end
else 
    Ex = Ex_short;
end

if iv
    if MSR
        Ey_short2 = reshape(Ey_short,[Nk_short,Nv]);
        Ey_short3 = Ey_short2.*Phase_short; 
        Ey = zeros(Nk,Nv);
        for ik = 1: Nk 
            Ey(ik,:) = ((p2(:,ik)).')*(Ey_short3.*Phase_long(ik,:));
        end
    else 
        Ey = Ey_short;
    end
    varargout{1} = Ey;
end

end