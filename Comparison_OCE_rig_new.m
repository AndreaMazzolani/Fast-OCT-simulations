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


saved_data = [path_str,'My_papers_code\Fast_OCT_simulations\git_hub\saved_data'];
str_sample = '/OCE_comparison_INCLUSION_Nk_433_Nscanx_409_Nscany_7_Nv_87500_Zmin_0_Zmax_500_zf_250.mat';

str_folder = [saved_data,'/',str_sample];

if not(isfolder(remote))&(exist(str_folder) == 2)
    DL = load(str_folder); if (isfield(DL,"DL")); DL = DL.DL; end
    unpackStruct(DL);
else

    %%%%%%%%%%%%%%%%%%%     Parameters of the simulation   %%%%%%%%%%%%%%%%%%%%%%%%
    % Scatterers region
    refind = 1;
    %     LX = 500e-6; maxZ =  500e-6; z_f = 300e-6/refind; str_image = 'ELL_C_scan'
    LX =  500e-6; maxZ =  500e-6; z_f = 300e-6/refind; str_image = 'ELL2_C_scan'
%     LX =  30e-6; maxZ =  100e-6; z_f = 10e-6/refind; str_image = 'ELL2_C_scan' %
    %     LX = 125e-6; maxZ = 1000e-6; z_f = 650e-6/refind; str_image = 'TEST_C_scan'
    %     LX = 750e-6; maxZ =  500e-6; z_f = 300e-6/refind; str_image = 'HARD'
    %     LX = 250e-6; maxZ =  500e-6; z_f = 250e-6/refind; str_image = 'INCLUSION'

    c = 2.997924580105029e+08;  % speed of light
    %     z_f = 300e-6/refind; 650e-6/refind;   % Axial coordinate of the focus of the Gaussian beam: 'z_f > 0'  -->  the mirror is placed at the left of the focal plane of the objective lens
    %     LX = 500e-6;125e-6;   % Maximum lateral position of the scatterers in x direction (meters)
    LY = 35e-6;     % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 000e-6;   % minimum of the scatterers axial position (meters)
    %     maxZ = 500e-6; 1000e-6;   % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 5; % Number of scatterers for ( 10 micron)^3
    Nv = round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    Nalpha_conf = 2;
    ALPHA_T = zeros(Nv,Nalpha_conf);
    % Lateral scanning p-arameters
    dscanx = 1.2e-6; Max_scanx = LX - 5e-6; Nscanx = ceil(2*Max_scanx/dscanx); % Number of A-scan for each B-scan
    dscany = 1.2e-6; Max_scany = LY - 5e-6; Nscany = 7; ceil(2*Max_scany/dscany); % Number of A-scan for each B-scan
%     dscanx = 3.2e-6; Max_scanx = 20e-6; LX - 5e-6; Nscanx = ceil(2*Max_scanx/dscanx); % Number of A-scan for each B-scan
%     dscany = 3.2e-6; Max_scany = 10e-6; LY - 5e-6; Nscany = 7; ceil(2*Max_scany/dscany); % Number of A-scan for each B-scan
    if (Nscanx == 1); x_scan_vec = 0e-6; else;  x_scan_vec = linspace(-Max_scanx,Max_scanx,Nscanx); end
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

    a_gr = 3e-3; a_smp = 0.5e-3;
    SIM = double(imread([path_str,'My_papers_code',slh,'Fast_OCT_simulations',slh,'saved_data',slh,str_image,'.png']));
    %     SIM = double(imread([path_str,'My_papers_code',slh,'Fast_OCT_simulations',slh,'saved_data',slh,'ELL_C_scan.png']));
    SIM = SIM(:,:,1); SIM = SIM./max(SIM(:)); if sum(isnan(SIM(:))) >0; error('The image file was not loaded properly!'); end
    SIM_ground = SIM.*0+a_gr; SIM = SIM_ground - SIM.*(a_gr-a_smp);
    Nx_sim = size(SIM,2); x_sim = linspace(-LX,LX,Nx_sim)'; dxs = x_sim(2)-x_sim(1);
    Nz_sim = size(SIM,1); z_sim = linspace(minZ,maxZ,Nz_sim)'; dzs = z_sim(2)-z_sim(1);
    ALPHA_SIM_tmp = v_U(:,1).*0;

    for ix = 1:Nx_sim
        fprintf('ix = %d / %d\n', ix,Nx_sim);
        x0 = x_sim(ix);
        for iz = 1:Nz_sim
            z0 = z_sim(iz); indxz = find( abs(z-z0) <= dzs/2 & abs(x-x0) <= dxs/2);
            ALPHA_SIM_tmp(indxz) = SIM(iz,ix);
        end
    end
    ALPHA_T(:,2) = ALPHA_SIM_tmp;     grid_size_alpha = [5,5,5]*1e-6;


    if 0;
        Nalpha_layers = 3;
        ALPHA = ones(Nalpha_conf, Nalpha_layers);% + (rand(Nalpha_conf, Nalpha_layers)-0.5).*dalpha;
        ALPHA(1,:) = [0,0,0]*1e-3; ALPHA(2,:) = [-1,1,3]*1e-3;
        Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
        borders = zeros(1,Nbord+2);        % Borders of the layers of strains
        borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
        borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
        borders(2:Nbord+1)  = minZ +  [200,400]*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
        BORD = ones(Nalpha_conf, Nalpha_layers + 1).*borders;


        for ib = 1:Nalpha_layers
            ibz = find((z >= borders(ib))&(z <= borders(ib+1)));
            ALPHA_T(ibz,2) = ALPHA(2,ib);
        end

        global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1;
        [alpha_output,IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor16(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,grid_size_alpha,ALPHA_T,x_scan_vec,y_scan_vec);

        global T_time LB NB ak_glob; T_time = 0; LB = 1; NB = 1; %ak_glob = [];
        [alpha_output,IT_CC2, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    end
    dx_sc = grid_size_alpha(1); dy_sc = grid_size_alpha(2); dz_sc = grid_size_alpha(3);
    v_L_TOT = load_scat_disp_3D(ALPHA_T,v_U,dx_sc,dy_sc,dz_sc);
    rho_T = rho(:).*ones(Nv,1,Nalpha_conf);

    global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1; %ak_glob = [];
    [IT_CC3, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor_gen(v_L_TOT,rho_T,f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);

    if 0;
        v_L = v_L_TOT(:,:,2);
        T_time = 0; LB = 1; NB = 1;
        [IU_CC_rig, IU_CC_AC_rig,IU_DC_sc_rig,IU_DC_ref_rig,akU_sc_rig,ak_ref_rig] = I_DW_Cscan(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);
        s_rig1 = seconds(T_time); s_rig.Format = 'hh:mm:ss'; s_rig;

        T_time = 0; LB = 1; NB = 1;
        [IL_CC_rig, IL_CC_AC_rig,IL_DC_sc_rig,IL_DC_ref_rig,akL_sc_rig,ak_ref_rig] = I_DW_Cscan(v_L,rho,f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);
        s_rig2 = seconds(T_time); s_rig.Format = 'hh:mm:ss'; s_rig;
        IT_CC_rig(:,:,:,1) = IU_CC_rig; IT_CC_rig(:,:,:,2) = IL_CC_rig;

        T_rig = s_rig1 + s_rig2;
        err_n(squeeze(IT_CC_rig(:,:,:,1)),squeeze(IT_CC(:,:,:,1)))
        comp(squeeze(IT_CC_rig(:,1,1,2)),squeeze(IT_CC(:,1,1,2)))
        err_n(IT_CC(:),IT_CC_rig(:))
    end

    % IU_CC = squeeze(IT_CC(:,:,:,1)); IL_CC = squeeze(IT_CC(:,:,:,2));
    % IU_CC2 = squeeze(IT_CC2(:,:,:,1)); IL_CC2 = squeeze(IT_CC2(:,:,:,2));
    IU_CC3 = squeeze(IT_CC3(:,:,:,1)); IL_CC3 = squeeze(IT_CC3(:,:,:,2));
    % IU_CC_rig = squeeze(IT_CC_rig(:,:,:,1)); IL_CC_rig = squeeze(IT_CC(:,:,:,2));

    f0 = mean(f_vec); f_max = max(f_vec);
    cut_off = 1e-3; Wk = 2*pi*(f_max-f0)/sqrt(log(1/cut_off)); Sk =  1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec(:)-f0)./Wk).^2);
    Df = f_vec(Nk)-f_vec(1); dz_tmp = 1e-6/sqrt(2); Nz = ceil(DZ/dz_tmp);
    [A_tmp,zb2] = gen_fft2(-1,f_vec(:),IU_CC3(:,1),2*[minZ,maxZ],Nz);  zbar = zb2/2; %zb2 = zbb2 - zbb2(Nz0); zb2 = zb2(:); zbar = zb2./2;
    % BU_CC = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2); BL_CC = gen_fft2(-1,f_vec,IL_CC.*Sk,zb2);
    % BU_CC2 = gen_fft2(-1,f_vec,IU_CC2.*Sk,zb2); BL_CC2 = gen_fft2(-1,f_vec,IL_CC2.*Sk,zb2);
    BU_CC3 = gen_fft2(-1,f_vec,IU_CC3.*Sk,zb2); BL_CC3 = gen_fft2(-1,f_vec,IL_CC3.*Sk,zb2);
    % BU_CC_rig = gen_fft2(-1,f_vec,IU_CC_rig.*Sk,zb2); BL_CC_rig = gen_fft2(-1,f_vec,IL_CC_rig.*Sk,zb2);
    str1 = [saved_data,'/OCE_comparison_',str_image,'_Nk_',num2str(Nk),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'.mat'];

    fit_range = 100e-6; res = 3.5e-6; Nx_lat = 4;  Ny_lat = 2;
    DL.x_scan_vec = x_scan_vec; DL.y_scan_vec = y_scan_vec; DL.Nscanx = Nscanx; DL.Nscany = Nscany;
    DL.Nx_lat = Nx_lat; DL.Ny_lat = Ny_lat; DL.a_gr = a_gr; DL.a_smp = a_smp;

    [str_xz, z_csr, Disp_xz,Err_xz,vec_lat_x,vec_lat_y] = CSR_Tay_10_tot_comb(minZ,maxZ,f_vec ,IU_CC3, IL_CC3,fit_range,res); Np_z = length(z_csr);
    DL.str_xz = str_xz; DL.z_csr = z_csr; DL.Disp_xz = Disp_xz; DL.Err_xz = Err_xz; DL.vec_lat_x = vec_lat_x; DL.vec_lat_y = vec_lat_y; save(str1,'DL');

    [str_xz2, z_csr2, Disp_xz2,Err_xz2,vec_lat_x2,vec_lat_y2] = CSR_Tay_10_tot_comb(minZ,maxZ,f_vec ,IU_CC3, IL_CC3,50e-6,res); Np_z2 = length(z_csr2);
    DL.str_xz2 = str_xz2; DL.z_csr2 = z_csr2; DL.Disp_xz2 = Disp_xz2; DL.Err_xz2 = Err_xz2; DL.vec_lat_x2 = vec_lat_x2; DL.vec_lat_y2 = vec_lat_y2; save(str1,'DL');

    [str_xz3, z_csr3, Disp_xz3,Err_xz3,vec_lat_x3,vec_lat_y3] = CSR_Tay_10_tot_comb(minZ,maxZ,f_vec ,IU_CC3, IL_CC3,fit_range,res,Nx_lat,Ny_lat,13,1); Np_z3 = length(z_csr3);
    DL.str_xz3 = str_xz3; DL.z_csr3 = z_csr3; DL.Disp_xz3 = Disp_xz3; DL.Err_xz3 = Err_xz3; DL.vec_lat_x3 = vec_lat_x3; DL.vec_lat_y3 = vec_lat_y3; save(str1,'DL');

    [str_xz3, z_csr3, Disp_xz3,Err_xz3,vec_lat_x3,vec_lat_y3] = CSR_Tay_10_tot_comb(minZ,maxZ,f_vec ,IU_CC3, IL_CC3,fit_range,res,Nx_lat,Ny_lat,13,1); Np_z3 = length(z_csr3);
    DL.str_xz3 = str_xz3; DL.z_csr3 = z_csr3; DL.Disp_xz3 = Disp_xz3; DL.Err_xz3 = Err_xz3; DL.vec_lat_x3 = vec_lat_x3; DL.vec_lat_y3 = vec_lat_y3; save(str1,'DL');

    iy0_rig = [4,5]; DL.iy0_rig = iy0_rig;
    BUL_tmp = BU_CC3.*conj(BL_CC3); BUL = permute(BUL_tmp,[2,3,1]); lambda0 = 1./f0;
    [str_UWA,z_UWA,unwrap_uwa_xyz] =  strain_UWA(lambda0,zbar,BUL(:,iy0_rig,:),1:Nscanx,1:2); %,1:Nscany);
    DL.str_UWA = str_UWA; DL.z_UWA = z_UWA;  save(str1,'DL');

    T_time = 0; LB = 1; NB = 1;
    [IU_CC_rig, IU_CC_AC_rig,IU_DC_sc_rig,IU_DC_ref_rig,akU_sc_rig,akU_ref_rig] = I_DW_Cscan(v_L_TOT(:,:,1),rho_T(:,:,1),f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec(iy0_rig));
    DL.IU_CC_rig = IU_CC_rig; DL.IU_CC3 = IU_CC3; DL.IL_CC3 = IL_CC3; save(str1,'DL');

    T_time = 0; LB = 1; NB = 1;
    [IL_CC_rig, IL_CC_AC_rig,IL_DC_sc_rig,IL_DC_ref_rig,akL_sc_rig,akL_ref_rig] = I_DW_Cscan(v_L_TOT(:,:,2),rho_T(:,:,2),f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec(iy0_rig));
    DL.IL_CC_rig = IL_CC_rig; save(str1,'DL');

    BU_CC_rig = gen_fft2(-1,f_vec,IU_CC_rig.*Sk,zb2); BL_CC_rig = gen_fft2(-1,f_vec,IL_CC_rig.*Sk,zb2);
    BUL_tmp = BU_CC_rig.*conj(BL_CC_rig); BUL_rig = permute(BUL_tmp,[2,3,1]);
    [str_UWA_rig,z_UWA_rig,unwrap_uwa_xyz] =  strain_UWA(lambda0,zbar,BUL_rig,1:Nscanx,iy0_rig);
    DL.str_UWA_rig = str_UWA_rig; DL.z_UWA_rig = z_UWA_rig; DL.iy0_rig = iy0_rig;
    DL.unwrap_uwa_xyz = unwrap_uwa_xyz; DL.BU_CC3 = BU_CC3; DL.BU_CC_rig = BU_CC_rig; DL.BL_CC3 = BL_CC3; DL.Nz = Nz; DL.zbar = zbar;  DL.x_scan_vec = x_scan_vec;

    [A_tmp,zb2_long] = gen_fft2(-1,f_vec(:),IU_CC3(:,1),2*[minZ,maxZ],Nz*11);  zbar_long = zb2_long/2; %zb2 = zbb2 - zbb2(Nz0); zb2 = zb2(:); zbar = zb2./2;
    BU_CC_rig_long = gen_fft2(-1,f_vec,IU_CC_rig.*Sk,zb2_long); BL_CC_rig_long = gen_fft2(-1,f_vec,IL_CC_rig.*Sk,zb2_long);
    BUL_tmp_long = BU_CC_rig_long.*conj(BL_CC_rig_long); BUL_rig_long = permute(BUL_tmp_long,[2,3,1]);
    [str_UWA_rig_long,z_UWA_rig_long,unwrap_uwa_xyz_long] =  strain_UWA(lambda0,zbar_long,BUL_rig_long,1:Nscanx,iy0_rig);
    DL.str_UWA_rig_long = str_UWA_rig_long; DL.z_UWA_rig_long = z_UWA_rig_long;  DL.zbar_long = zbar_long;
    DL.unwrap_uwa_xyz_long = unwrap_uwa_xyz_long; DL.BUL_rig_long = BUL_rig_long;
    save(str1,'DL');

    disp(str1)

end

if not(strcmp(path_str,remote))
    str_UWA2 = permute(str_UWA,[3,1,2]); str_UWA_rig2 = permute(str_UWA_rig,[3,1,2]); str_UWA_rig_long2 = permute(str_UWA_rig_long,[3,1,2]);

    Nscany_0 = length(iy0_rig);
    Nx_slit = min(6,Nscanx); Ny_slit = min(3,Nscany_0); Nz_slit = 3; Nz_slit_long = Nz_slit*11;
    str_xz_opt = strain_optimized(str_xz,Err_xz,Nx_slit,Ny_slit,Nz_slit,22);
    str_xz2_opt = strain_optimized(str_xz2,Err_xz2,Nx_slit,Nscany,Nz_slit,22);
    str_xz3_opt = strain_optimized(str_xz3,Err_xz3,Nx_slit,Ny_slit-2*Ny_lat,Nz_slit,22);
    str_uwa_opt = strain_optimized(str_UWA2,str_UWA2.*0+1,Nx_slit,Ny_slit,Nz_slit,22);
    str_UWA_rig_opt = strain_optimized(str_UWA_rig2,str_UWA_rig2.*0+1,Nx_slit,Ny_slit,Nz_slit,22);
    str_UWA_rig_long_opt = strain_optimized(str_UWA_rig_long2,str_UWA_rig_long2.*0+1,Nx_slit,Ny_slit,Nz_slit_long,22);

    figure(1); clf; NF = 13;  sgtitle('Strain retrieved (WLS)', 'FontSize',2*NF); da = 0.3e-3; a_min = -max(a_gr,a_smp) - da; a_max = -min(a_gr,a_smp) +2* da; vec_a = [a_min,a_max];
    subplot(2,1,1); imagesc_set(str_UWA_rig_opt,z_UWA_rig(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('rigorous'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_a); colormap('hot');  axis equal;
    subplot(2,1,2); imagesc_set(str_uwa_opt,z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('Taylor'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_a);  colormap('hot');   axis equal;

    res_str = str_uwa_opt-str_UWA_rig_opt; abs_err = abs(res_str);
    figure(12); clf;  NF = 13; ixb = ceil(Nscanx.*0.4); MB = max(abs(BU_CC3(:,:,iy0_rig(1))),[],'all'); vec_b = [-2,0]; sgtitle('Strain retrieved (WLS)', 'FontSize',2*NF);
    vec_res = [0,1]*1e-4; da = 0.3e-3; a_min = -max(a_gr,a_smp) - da; a_max = -min(a_gr,a_smp) +2* da; vec_a = [-4,1]*1e-3;
    subplot(2,2,1); imagesc_set(str_UWA_rig_opt,z_UWA_rig(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('rigorous'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_a); colormap('hot');  axis equal;
    subplot(2,2,2); imagesc_set(str_uwa_opt,z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('Taylor'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_a);  colormap('hot');   axis equal;
    subplot(2,2,3); imagesc_set(abs_err,z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('Residual image'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_res);  colormap('hot');   axis equal;
    subplot(2,2,4); semilogy(x_scan_vec(Nx_slit:end).*1e6,max(abs_err)); ylim([3e-8,3e-4]); title('Maximum error'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF);

    figure(3); clf; vec_a = [-1,0]; MU = max(abs(BU_CC_rig(:,:,1)),[],'all'); BU_app = log10(abs(BU_CC3(:,:,iy0_rig(1))./MU));
    BU_rig = log10(abs(BU_CC_rig(:,:,1)./MU)); BU_res =  log10(abs((BU_CC_rig(:,:,1)./MU)-(BU_CC3(:,:,iy0_rig(1))./MU)));
    subplot(2,1,1); imagesc_set((BU_rig),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_a);
    %     subplot(3,1,2); imagesc_set((BU_app),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_a);
    subplot(2,1,2); imagesc_set(BU_res,zbar*1e6,x_scan_vec*1e6); axis equal; caxis([-5,-3]); colormap('gray');

    err_threshold_x = x_scan_vec.*0+ log10(3e-4);
    figure(113); clf; vec_b = [-1.5,0]; vec_r1 = [-5,-3]; vec_r2 = [-5,-2]; vec_s = [-3.5,0]*1e-3;
    vec_rs1 = [-8,-4]; vec_rs2 = [-6,-2];  x_str = [x_scan_vec(3+Nx_slit),x_scan_vec(end)]*1e6;
    subplot(2,4,1); imagesc_set((BU_rig),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_b);
    subplot(2,4,2); imagesc_set((BU_app),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_b);
    subplot(2,4,3); imagesc_set((BU_res),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_r1);
    subplot(2,4,4); plot(x_scan_vec*1e6,log10(2*sum(abs(BU_CC3(:,:,iy0_rig(1))-BU_CC_rig(:,:,1)))./(sum(abs(BU_CC3(:,:,iy0_rig(1))))+sum(abs(BU_CC_rig(:,:,1))))),'k'); hold on; plot(x_scan_vec.*1e6,err_threshold_x,'r--'); ylim(vec_r2); xlim(x_str);
    %     subplot(2,4,4); plot(x_scan_vec*1e6,log10(max(abs(BU_CC3(:,:,iy0_rig(1))./MU-BU_CC_rig(:,:,1)./MU),[],1)),'k'); ylim(vec_r2); xlim(x_str);
    subplot(2,4,5); imagesc_set(str_UWA_rig_opt,z_UWA_rig(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('rigorous'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_s); colormap('gray');  axis equal;
    subplot(2,4,6); imagesc_set(str_uwa_opt,z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('app'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_s); colormap('gray');  axis equal;
    subplot(2,4,7); imagesc_set(log10(abs_err),z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('Residual image'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_rs1);     axis equal;
    subplot(2,4,8); plot(x_scan_vec(Nx_slit:end)*1e6,log10(2*sum(abs(res_str))./(sum(abs(str_uwa_opt))+sum(abs(str_UWA_rig_opt)))),'k');  hold on; plot(x_scan_vec.*1e6,err_threshold_x,'r--'); ylim(vec_rs2); xlim(x_str);
    %     subplot(2,4,8); plot(x_scan_vec(Nx_slit:end)*1e6,log10(max(res_str,[],1)),'k'); ylim(vec_rs2); xlim(x_str);

    figure(81); clf; imagesc_set((BU_rig),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_b);colormap('gray');
    figure(82); clf; imagesc_set((BU_app),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_b);colormap('gray');
    figure(83); clf; imagesc_set((BU_res),zbar*1e6,x_scan_vec*1e6); axis equal; caxis(vec_r1);colormap('turbo');
    figure(84); clf; plot(x_scan_vec*1e6,log10(2*sum(abs(BU_CC3(:,:,iy0_rig(1))-BU_CC_rig(:,:,1)))./(sum(abs(BU_CC3(:,:,iy0_rig(1))))+sum(abs(BU_CC_rig(:,:,1))))),'k'); hold on; plot(x_scan_vec.*1e6,err_threshold_x,'r--');colormap('gray'); ylim(vec_r2); xlim(x_str);
    figure(85); clf; imagesc_set(str_UWA_rig_opt,z_UWA_rig(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('rigorous'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_s); colormap('hot');  axis equal;
    figure(86); clf; imagesc_set(str_uwa_opt,z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('app'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_s); colormap('hot');  axis equal;
    figure(87); clf; imagesc_set(log10(abs_err),z_UWA(Nz_slit:end)*1e6,x_scan_vec(Nx_slit:end)*1e6); title('Residual image'); ylabel('z(\mu m)'); xlabel('x(\mu m)'); set(gca, 'FontSize',NF); caxis(vec_rs1); colormap('hot');    axis equal;
    figure(88); clf; plot(x_scan_vec(Nx_slit:end)*1e6,log10(2*sum(abs(res_str))./(sum(abs(str_uwa_opt))+sum(abs(str_UWA_rig_opt)))),'k');  hold on; plot(x_scan_vec.*1e6,err_threshold_x,'r--'); ylim(vec_rs2); xlim(x_str);

end


