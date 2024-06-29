clear all; clc;

slh = '\'; if isunix; slh = '/'; end;

desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus
laptop = 'C:/Users/andre/Documents/MATLAB/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus

if isfolder(desktop)
    path_str = desktop;
elseif isfolder(laptop)
    path_str = laptop;
elseif isfolder(remote)
    path_str = remote;
end

rmpath(genpath(path_str));
addpath(genpath([path_str,'My_papers_code\Fast_OCT_simulations\git_hub\']));



refind = 1;
LX = 50e-6; maxZ =  300e-6; z_f = 50e-6/refind; str_image = 'ELL2_C_scan'

c = 2.997924580105029e+08;  % speed of light
LY = 35e-6;     % Maximum lateral position of the scatterers in y direction (meters)
minZ = 000e-6;   % minimum of the scatterers axial position (meters)
DZ = maxZ - minZ; % Width of the axial interval
nn3_dens = 5; % Number of scatterers for ( 10 micron)^3
Nv = round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

% Lateral scanning p-arameters
dscanx = 3.2e-6; Max_scanx = LX - 5e-6; Nscanx = ceil(2*Max_scanx/dscanx); % Number of A-scan for each B-scan
dscany = 3.2e-6; Max_scany = LY - 5e-6;  ceil(2*Max_scany/dscany); % Number of A-scan for each B-scan
if (Nscanx == 1); x_scan_vec = 0e-6; else; Nscany = 1; x_scan_vec = linspace(-Max_scanx,Max_scanx,Nscanx); end
if (Nscany == 1); y_scan_vec = 0e-6; else;  y_scan_vec = linspace(-Max_scany,Max_scany,Nscany); end

W_0 = 4.6e-6; NA = 0.0972; foc1 = 25*1e-3; foc2 = 36*1e-3; b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)

%%%%%%%%%%%%%%    Unloaded scatterer region setting    %%%%%%%%%%%%%%%%%%%%%%
rhoMax = 3e-4; rho = ones(1,Nv).*rhoMax;
v_U = zeros(Nv,3); % Array of scatterers initialization: 2D array of the scatterers in the unloaded case. This is an array of dimension (Nv,3), where each row consists of 3 numbers that are the (x,y,z) coordinates of each scatterer
v_U(:,1) = -LX + rand(Nv,1)*2*LX; % (x coordinates of all scatterers)
v_U(:,2) = -LY + rand(Nv,1)*2*LY; % (y coordinates of all scatterers)
v_U(:,3) = sort(minZ + rand(Nv,1)*DZ);  % (z coordinates of all scatterers)
x = v_U(:,1); y = v_U(:,2); z = v_U(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%     Spatial frequency setting  'f = 1/lambda'    %%%%%%%
lam_min = 1170.5e-9; lam_max = 1407.8e-9; f_min = 1./lam_max; f_max = 1./lam_min;
f0 = ( f_max + f_min)/2; Df = f_max - f_min; Nk_nyq = 4*max(abs(v_U(:,3).*refind))*Df;
Nk =  ceil(Nk_nyq*1.5); % Minimum number of frequencies able to simulate the deepest axial position (coming from the Nyquist frequency and rates relation)% Number of wavenumbers, they are set as 1.5 times 'Nk_nyq' ( A safe number greater than 'Nk_nyq')
f_vec = linspace(f_min,f_max,Nk); f_vec_t = f_vec.*c;

Nalpha_conf = 10; Nalpha_layers = 2; 
ALPHA =  zeros(Nalpha_conf, Nalpha_layers);
alpha_vec = linspace(0,5,Nalpha_conf)*1e-3;
ALPHA(:,1) = alpha_vec(:); % Array of strains for "version 15"
ALPHA(:,2) = 1.5*alpha_vec(:); % Array of strains for "version 15"
  
ALPHA_T = zeros(Nv,Nalpha_conf); % Array of strains for "version 16"


Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
borders = zeros(1,Nbord+2);        % Borders of the layers of strains
borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
borders(2:Nbord+1)  = minZ + [150]*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
BORD = ones(Nalpha_conf, Nalpha_layers + 1).*borders;

for ia = 1:Nalpha_conf
    for ic = 1:Nalpha_layers; 
        IND_tmp = find( (v_U(:,3)>BORD(ia,ic)) & ( v_U(:,3)<=BORD(ia,ic+1)) );
        ALPHA_T(IND_tmp,ia) = ALPHA(ia,ic);
    end
end

grid_size_alpha = [5,5,5]*1e-6;
dx_sc = grid_size_alpha(1); dy_sc = grid_size_alpha(2); dz_sc = grid_size_alpha(3);
v_L_TOT = load_scat_disp_3D(ALPHA_T,v_U,dx_sc,dy_sc,dz_sc);
rho_T = rho(:).*ones(Nv,1,Nalpha_conf);

global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1; 
[alpha_output,IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
time_15 = T_time;

global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1;
[alpha_output2,IT_CC2, IT_CC_AC2,IT_DC_sc2,IT_DC_ref2,akT_sc2,akT_ref2,ERROR_mat2,ak_ref2] = Series_Loaded_DW_Cscan_uniform_Taylor16(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,grid_size_alpha,ALPHA_T,x_scan_vec,y_scan_vec);
time_16 = T_time;

global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1;  
[IT_CC3, IT_CC_AC3,IT_DC_sc3,IT_DC_ref3,akT_sc3,akT_ref3,ERROR_mat3,ak_ref3] = Series_Loaded_DW_Cscan_uniform_Taylor_gen(v_L_TOT,rho_T,f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);
time_gen = T_time;


fprintf('Computational time for version:\n 15 : %f secs ,\n 16 : %f secs ,\n gen : %f secs \n',time_15,time_16,time_gen)
fprintf('Rel. error: 15-16: %f \n 15-gen: %f \n gen-16: %f\n',err_n(IT_CC,IT_CC2),err_n(IT_CC,IT_CC3),err_n(IT_CC2,IT_CC3))

figure(1); clf; hold on; plot(IT_CC(:,14,1,1));  plot(IT_CC2(:,14,1,1)); plot(IT_CC3(:,14,1,1));
figure(2); clf; hold on; plot(IT_CC(:,14,1,3));  plot(IT_CC2(:,14,1,3)); plot(IT_CC3(:,14,1,3));

cut_off = 1e-3; Wk = 2*pi*(f_max-f0)/sqrt(log(1/cut_off)); Sk =  1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec(:)-f0)./Wk).^2);
Df = f_vec(Nk)-f_vec(1); dz_tmp = 1e-6/sqrt(2); Nz = ceil(DZ/dz_tmp);

[BT_CC,zb2] = gen_fft2(-1,f_vec,IT_CC.*Sk,2*[minZ,maxZ],Nz);  zbar = zb2/2;  BU_CC = BT_CC(:,:,:,1);
BT_CC2 = gen_fft2(-1,f_vec,IT_CC2.*Sk,zb2);  BU_CC2 = BT_CC2(:,:,:,1);
BT_CC3 = gen_fft2(-1,f_vec,IT_CC3.*Sk,zb2);   BU_CC3 = BT_CC3(:,:,:,1);

ia0 = 6;
IU_CC = squeeze(IT_CC3); IL_CC = squeeze(IT_CC(:,:,:,ia0));
IU_CC2 = squeeze(IT_CC2); IL_CC2 = squeeze(IT_CC2(:,:,:,ia0));
IU_CC3 = squeeze(IT_CC3); IL_CC3 = squeeze(IT_CC3(:,:,:,ia0));

% str1 = [DW_signal_folder,'/OCE_comparison_',str_image,'_Nk_',num2str(Nk),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'.mat'];

fit_range = 100e-6; res = 3.5e-6; Nx_lat = 4;  Ny_lat = 2; lambda0 = 1./f0;

BUL_tmp = BU_CC .*conj(BT_CC ); BUL  = permute(BUL_tmp,[2,3,1,4]); 
BUL_tmp = BU_CC2.*conj(BT_CC2); BUL2 = permute(BUL_tmp,[2,3,1,4]); 
BUL_tmp = BU_CC3.*conj(BT_CC3); BUL3 = permute(BUL_tmp,[2,3,1,4]); 

clear str_UWA str_UWA2 str_UWA3
for ia = 2:Nalpha_conf
    ia
    [str_UWA(:,:,:,ia),z_UWA,unwrap_uwa_xyz] =  strain_UWA(lambda0,zbar,BUL(:,:,:,ia),1:Nscanx,1:Nscany);
    [str_UWA2(:,:,:,ia),z_UWA,unwrap_uwa_xyz] =  strain_UWA(lambda0,zbar,BUL2(:,:,:,ia),1:Nscanx,1:Nscany);
    [str_UWA3(:,:,:,ia),z_UWA,unwrap_uwa_xyz] =  strain_UWA(lambda0,zbar,BUL3(:,:,:,ia),1:Nscanx,1:Nscany);
end

sA = permute(str_UWA,[3,1,4,2]); sA2 = permute(str_UWA2,[3,1,4,2]); sA3 = permute(str_UWA3,[3,1,4,2]);

figure(3); clf; ia0 = 4; ix0 = 23;
plot(z_UWA*1e6,sA(:,ix0,ia0),z_UWA*1e6,sA2(:,ix0,ia0),z_UWA*1e6,sA3(:,ix0,ia0),'--'); legend('15','16','gen');



figure(4); clf;  veca = -flip([alpha_vec(ia0), alpha_vec(ia0)*1.5]);
subplot(2,2,1); imagesc_set(sA(:,:,ia0)); caxis(veca);
subplot(2,2,2); imagesc_set(sA2(:,:,ia0)); caxis(veca);
subplot(2,2,3); imagesc_set(sA3(:,:,ia0)); caxis(veca);
% subplot(2,2,4); imagesc_set(sAr(:,:,ia0)); sAr needs to be made

