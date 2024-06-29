% Simulation of a Large OCT C-scan with the Taylor expansion function "15"

clear all; clc;   % Inizialization

TEST = 1; slh = '\'; if isunix; slh = '/'; end
laptop = 'C:\Users\andre\Documents\MATLAB\Matlab_drive_scripts\'; % Desktop personal folder
desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus

if isfolder(laptop)
    path_str = laptop;
elseif isfolder(desktop)
    path_str = desktop;
elseif isfolder(remote)
    path_str = remote;
end

addpath(genpath([path_str,'handle_functions']));
addpath(genpath([path_str,'Illumination_functions']));
addpath(genpath([path_str,'saved_heavy_files']));
addpath(genpath([path_str,'OCT_models']));

% Set the Matlab current folder where the folder --> "Simulated_OCT_DW_UL_Bscan" was copied
% cd([path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan'])

tic;

% Folders to add to the path
if isfolder(remote)
    DW_signal_folder = [path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan/DW_seriesL_uniform_signals'];
else
    DW_signal_folder = [path_str,'saved_heavy_files/fast_code_paper'];
end
% DW_signal_folder = [path_str,'saved_heavy_files'];
% DW_signal_folder = 'saved_data';

str_folder = [DW_signal_folder,'/' ...%dfgag'];%
  'Large_OCT_Cscan_new_Nk_865_W0_9_Nscanx_200_Nscany_200_Nv_312500_Zmin_0_Zmax_1e+03_zf_290_alpha_0_0_NL_1_NA_1.mat'];
%   'Large_OCT_Cscan_new_Nk_865_W0_9_Nscanx_200_Nscany_3_Nv_312500_Zmin_0_Zmax_1e+03_zf_290_alpha_0_0_NL_1_NA_1.mat'];
if not(isfolder(remote)) & (exist(str_folder) == 2)  
    DL = load(str_folder); if (isfield(DL,"DL")); DL = DL.DL; end; unpackStruct(DL);
else
    answer = input('There are not loaded files, do you still want to proceed?','s');
    while not((answer == "n")|(answer == "y"))
        answer = input(['Answer y or n:\n'],'s');
    end
    if answer == "n"
        error(["Check the file name!"]);
    end

    c = 2.997924580105029e+08;  % speed of light

    %%%%%%%%%%%%%%%%%%%     Parameters of the simulation   %%%%%%%%%%%%%%%%%%%%%%%%
    % Scatterers region
    refind = 1;1.42; 
    z_f = 290e-6/refind; % Axial coordinate of the focus of the Gaussian beam: 'z_f > 0'  -->  the mirror is placed at the left of the focal plane of the objective lens
    z_f2 = 750e-6/refind;
    LX = 125e-6;     % Maximum lateral position of the scatterers in x direction (meters)
    LY = 125e-6;     % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 000e-6;   % minimum of the scatterers axial position (meters)
    maxZ = 1000e-6;  % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 5; % Number of scatterers for ( 10 micron)^3
    Nv =  round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    dalpha = 11e-4;
    alpha_min = 0e-3;
    alpha_max = 0e-3;
    NA_tot = 1;
    Nalpha_conf = 1;NA_tot*2+1;
    Nalpha_layers = 1; % Number of uniform layers for each combination
    alpha_VEC = 0; linspace(alpha_min + 0.*dalpha/2,alpha_max - 0.*dalpha/2,Nalpha_conf).';% alpha_min + rand(Nalpha_conf,Nalpha_layers) * (alpha_max-alpha_min);
    %   ALPHA = ones(Nalpha_conf, Nalpha_layers).*alpha_VEC;% + (rand(Nalpha_conf, Nalpha_layers)-0.5).*dalpha;
    ALPHA = zeros(Nalpha_conf, Nalpha_layers);

    %     ALPHA(2:NA_tot+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(2:NA_tot+1,2) = ALPHA(2:NA_tot+1,1);
    %     ALPHA(NA_tot+2:NA_tot*2+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(NA_tot+2:NA_tot*2+1,2) = 1e-3.*10;
    Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
    borders = zeros(1,Nbord+2);        % Borders of the layers of strains
    borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
    borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
    borders(2:Nbord+1)  = minZ +  [1000]*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
    BORD = ones(Nalpha_conf, Nalpha_layers + 1).*borders;


    % Control statements to check that the array of borders has been set correctly.
    if (min(min(diff(BORD,[],2))) < 0);
        error(["Borders must be an increasing array! Borders coordinates originate at the focal plane. "])
    end
    if sum((BORD(1,:) < minZ)) > 1;
        error(["First border must be less or equal than the shallowest scatterers "]);
    end

    % Reflectivity coefficients of the scattererss
    rhoMax = 3e-4; % Maximum of the reflectivity profile (This is a relative parameter, considering the mirror reflectivity = 1!!);

    % Lateral scanning p-arameters
    dscan = 1.2e-6;
    Max_scan = LX -5e-6;- 30e-6; % Maximum lateral scanning in positive and negative x-direction (meters)
    Nscanx =   ceil(2*Max_scan/dscan); % Number of A-scan for each B-scan
    if (Nscanx == 1); x_scan_vec = 0e-6; y_scan_vec = -0e-6; Nscany = Nscanx; else;  x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx); y_scan_vec = 0; Nscany = length(y_scan_vec); end;
    % Nscany = 1; y_scan_vec = 0; Nscanx0 = ceil(Nscanx/2)+1; Nscanx = 2*Nscanx0-1; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx);
    % Uniform refractive index of the medium
    y_scan_vec =  x_scan_vec;[-1:1]*1e-6; Nscany = length(y_scan_vec);

    % The following physical parameters have been set as the OCT at UCL
    W_0 = 9e-6;7e-6;11.6e-6; % Waist radius of the Gaussian beam at the aperture plane between the collimator and the objective lens
    NA = 0.0972;   % Numerical aperture (Ra = NA*foc2 --> NA = Ra/foc2)
    %   Ra = 3.5e-3; % Radius of the aperture
    %   NA = Ra/foc2;
    lam_min = 1170.5e-9;  % minmum wavelength (meters)
    lam_max = 1407.8e-9;  % maximum wavelength (meters)

    foc1 = 25*1e-3; % focal length of the collimator (meters)
    foc2 = 36*1e-3; % focal length of the objective lens (meters)
    %     Ra = 3.5e-3;

 

    % Dependent parameters
%     rho = rand(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    rho = ones(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%    Unloaded scatterer region setting    %%%%%%%%%%%%%%%%%%%%%%
    v_U = zeros(Nv,3); % Array of scatterers initialization: 2D array of the scatterers in the unloaded case. This is an array of dimension (Nv,3), where each row consists of 3 numbers that are the (x,y,z) coordinates of each scatterer
    v_U(:,1) =   -LX + rand(Nv,1)*2*LX; % (x coordinates of all scatterers)
    v_U(:,2) =   -LY + rand(Nv,1)*2*LY; % (y coordinates of all scatterers)
    v_U(:,3) =   sort(minZ + rand(Nv,1)*DZ);  % (z coordinates of all scatterers)
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SIM = double(imread([path_str,'My_papers_code',slh,'Fast_OCT_simulations',slh,'saved_data',slh,'TEST_C_scan.png']));
    SIM = SIM./max(SIM(:)); SIM = SIM(:,:,1); 
    Nx_sim = size(SIM,2); x_sim = linspace(-LX,LX,Nx_sim)'; dxs = x_sim(2)-x_sim(1);
    Nz_sim = size(SIM,1); z_sim = linspace(0,maxZ,Nz_sim)'; dzs = z_sim(2)-z_sim(1);
    rho_SIM_tmp = rho.*0;  x = v_U(:,1); y = v_U(:,2); z = v_U(:,3);

    for ix = 1:Nx_sim
        fprintf('ix = %d / %d\n', ix,Nx_sim)
        x0 = x_sim(ix);  
        for iz = 1:Nz_sim
            z0 = z_sim(iz); indxz = find( abs(z-z0) <= dzs/2 & abs(x-x0) <= dxs/2);
            rho_SIM_tmp(indxz) = SIM(iz,ix);
        end
    end
% save([path_str,'My_papers_code\Fast_OCT_simulations\saved_data\rho_SIM_tmp.mat'],'rho_SIM_tmp');
rho_SIM = rho_SIM_tmp*rhoMax*10 + rho;
rho_SIM2 = rho_SIM_tmp*rhoMax*3 + rho;

    %%%%%%%%%   To Test the model   %%%%%%%%%
    % xg = linspace(min(x), max(x), 301);
    % zg = linspace(min(z), max(z), 1201);
    % [Xg, Zg] = meshgrid(xg, zg);
    % rhog = griddata(x, z, rho_SIM, Xg, Zg);
    % 
    % rhog(isnan(rhog)) = 0;
    % clf; imagesc_set(rhog)

%     ind_sp = find(x.^2+(z-z_f).^2 <= (10e-6).^2); rho(ind_sp) = rhoMax*10;
    %%%%%%%%%%%%%%     Spatial frequency setting  'f = 1/lambda'    %%%%%%%
    f_min = 1./lam_max;  % minimum frequency of the simulated spectrum
    f_max = 1./lam_min;  % maximum frequency of the simulated spectrum
    f0 = ( f_max + f_min)/2; % central frequency of the simulated spectrum
    Df = f_max - f_min; % frequencies spectrum width
    Nk_nyq = 4*max(abs(v_U(:,3).*refind))*Df; % Minimum number of frequencies able to simulate the deepest axial position (coming from the Nyquist frequency and rates relation)

    Nk =  ceil(Nk_nyq*1.5); % Number of wavenumbers, they are set as 1.5 times 'Nk_nyq' ( A safe number greater than 'Nk_nyq')
    f_vec = linspace(f_min,f_max,Nk);  % Spatial-frequencies array
    f_vec_t = f_vec.*c; % Time-frequencies_array
    k_vec = 2*pi*f_vec; % Arrays of wavenumbers;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scattered field and modal coefficients calculation (First Born approximation)
    str_a = [num2str(alpha_min),'_',num2str(alpha_max),'_NL_',num2str(Nalpha_layers),'_NA_',num2str(Nalpha_conf)];
    str = [DW_signal_folder,'/Large_OCT_Cscan_new_Nk_',num2str(Nk),'_W0_',num2str(W_0*1e6),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'_alpha_',str_a,'.mat'];

    global T_time LB NB;
    %x_scan_vec2 = linspace(x_scan_vec(1),x_scan_vec(end),301);
    T_time = 0; LB = 1; NB = 1; %ix_sub = 100:150; 
    [alpha_output,IT_0, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat]  = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);  DL.IT_0 = IT_0;      save(str,'DL');
    [alpha_output,IT_1, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat]  = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho_SIM_tmp*rhoMax,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec); DL.IT_1 = IT_1;   save(str,'DL'); 
    [alpha_output,IT_z0, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f2,BORD,ALPHA,x_scan_vec,y_scan_vec); DL.IT_z0 = IT_z0;   save(str,'DL');
    [alpha_output,IT_z1, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho_SIM_tmp(:)*rhoMax,f_vec_t,refind,W_0/b12,NA,z_f2,BORD,ALPHA,x_scan_vec,y_scan_vec); DL.IT_z1 = IT_z1;   save(str,'DL');
%     [alpha_output,IT_z2, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat] = Series_Loaded_DW_Cscan_uniform_Taylor_gen(v_U,rho_SIM_tmp(:)*rhoMax,f_vec_t,refind,W_0/b12,NA,z_f2,x_scan_vec(ix_sub),y_scan_vec);
    s = seconds(T_time); s.Format = 'hh:mm:ss'; simulation_time = s; fprintf(['\nMax error s = ',num2str(max(ERROR_mat)), ', Computational time = ']); 

    % subset of data for the test
    ik_sub = 1:Nk; ix_sub = ceil(Nscanx*0.5); iy_sub = ceil(Nscany*0.5); 
    T_time = 0; LB = 1; NB = 1; global ak_glob; ak_glob = squeeze(akT_ref(ik_sub,1,1));
    [IT_CC_rig, IT_CC_AC_rig,IT_DC_sc_rig,IT_DC_ref_rig,ak_sc_rig,ak_ref_rig] = I_DW_Cscan(v_U,rho,f_vec_t(ik_sub),refind,W_0/b12,NA,z_f,x_scan_vec(ix_sub),y_scan_vec(iy_sub));
    t_rig = T_time; s_rig = seconds(T_time); s_rig.Format = 'hh:mm:ss'; s_rig
    s_rig_single_A_scan = s_rig;  DL.s_rig_single_A_scan = s_rig_single_A_scan;  save(str,'DL');

    
    DL.f_vec = f_vec; DL.Nscanx = Nscanx;  DL.Nscany = Nscany;  DL.Nk = Nk; DL.refind = refind; 
    DL.z_f2 = z_f2; DL.z_f = z_f; DL.maxZ = maxZ; DL.x_scan_vec = x_scan_vec; DL.y_scan_vec = y_scan_vec;
    DL.simulation_time = simulation_time;

    save(str,'DL'); disp(str);

end

if not(isfolder(remote))
    Estimated_rigorous_time = Nscanx*Nscany*s_rig_single_A_scan;
    disp(['Estimated time for the rigorous calculation: ']); Estimated_rigorous_time.Format = 'dd:hh:mm:ss'

    f0 = mean(f_vec); f_max = max(f_vec); f_min = min(f_vec); Df = f_max-f_min;
    cut_off =  2e-1; Wk = 2*pi*(max(f_vec)-f0)/sqrt(log(1/cut_off)); Sk =  hanning(Nk); ;1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec(:)-f0)./Wk).^2); hanning(Nk); 
    Nyq_z = Nk/Df/4; zmin = 0;-Nyq_z; zmax = maxZ*refind + 50e-6; Nyq_z; Nz = 1101; Nx0 = ceil(Nscanx/2);
    [A_tmp2,zbb2] = gen_fft2(-1,f_vec(:)-f0,Sk,2*[zmin,zmax],Nz); Nz0 = min(find(zbb2 >= 0));   zb2 = zbb2 -zbb2(Nz0); zbar = zb2/2;

    IU_CC = squeeze(IT_0);  BU_CC2 = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2);
    IU_CC = squeeze(IT_1);  BU_CC1 = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2);
    IU_CC = squeeze(IT_z0);  BU_CCz2 = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2);
    IU_CC = squeeze(IT_z1);  BU_CCz1 = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2);

    BU_CC = BU_CC1*5+BU_CC2; BU_CCz = BU_CCz1*5+BU_CCz2; coef_db = 0.2;
    MA = mean( (abs(BU_CC(:)))).*coef_db;   MAz = mean( (abs(BU_CCz(:)))).*coef_db; 
    IZ_min = min(find(zbar > 0e-6)); IZ_max = max(find(zbar < 1000e-6)); IZ = IZ_min:IZ_max; Nz2 = length(IZ);
    
    figure (1); clf; ca1 = 0; ca2 = 66; %print(['fig1'],'-dpdf','-bestfit')
    subplot(1,2,1); imagesc_set(abs(BU_CC(IZ,:,3)./MA),zbar(IZ),x_scan_vec); axis equal; caxis([ca1,ca2]);   colormap('gray')
    subplot(1,2,2); imagesc_set(abs(BU_CCz(IZ,:,3)./MAz),zbar(IZ),x_scan_vec); axis equal; caxis([ca1,ca2]);  colormap('gray')

    figure (11); clf; ca1 = 0; ca2 = 30; corr = 2; %print(['fig11'],'-dpdf','-bestfit')
    subplot(1,2,1); imagesc_set(20*log10(abs(BU_CC(IZ,:,3))./MA),zbar(IZ),x_scan_vec); axis equal;   colormap('gray'); caxis([ca1,ca2]); 
    subplot(1,2,2); imagesc_set(20*log10(abs(BU_CCz(IZ,:,3)./MAz)),zbar(IZ),x_scan_vec); axis equal;   colormap('gray'); caxis([ca1,ca2]);

    Nv_tot = Nz2*Nscanx*Nscany; v_oct = zeros(Nv_tot,3); B_oct = zeros(Nv_tot,1); B_octz = zeros(Nv_tot,1);
    
    
    for ix = 1:Nscanx; %1:Nscanx 
        ix
        for iy = 1:Nscany
            it1 = 1 + Nz2*((iy-1) + (ix-1)*Nscany); it2 = it1 + Nz2-1;
            v_oct(it1:it2,:) = [x_scan_vec(ix)*ones(Nz2,1),y_scan_vec(iy)*ones(Nz2,1),flip(zbar(IZ))']*1e6;
            B_oct(it1:it2)  = abs(BU_CC(IZ,ix,iy))./MA;
            B_octz(it1:it2) = abs(BU_CCz(IZ,ix,iy))./MAz;
        end
    end

%     ptCloud = pointCloud(v_oct); ptCloud.Intensity = B_oct;
%     ptCloudz = pointCloud(v_oct); ptCloud2.Intensity = B_octz;

%     figure(2); clf; NF = 11;
%     subplot(1,2,1); pcshow(ptCloud); set(gcf,'color','w'); colorbar; set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]); colormap('gray');
%     subplot(1,2,2); pcshow(ptCloudz); set(gcf,'color','w'); colorbar;  set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);

    figure(2); clf; NF = 11; sgtitle('|C-scan|','FontSize',NF*2); print(['fig2'],'-dpdf','-bestfit')
    subplot(1,2,1); pcshow(v_oct,log10(B_oct)); set(gcf,'color','w'); set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]); colormap('gray');  xlabel('x (\mu m)');  ylabel('y (\mu m)');  zlabel('z (\mu m)');  title('$z_f = 290 \mu m$'); % caxis([ca1,ca2]); 
    subplot(1,2,2); pcshow(v_oct,(log10(B_octz))); set(gcf,'color','w'); colorbar;  set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]); colormap('gray');  xlabel('x (\mu m)');  ylabel('y (\mu m)');  zlabel('z (\mu m)'); title('$z_f = 750 \mu m$'); % caxis([ca1,ca2]);

    figure(12); clf; NF = 11; sgtitle('|C-scan|','FontSize',NF*2); print(['fig12'],'-dpdf','-bestfit')
    subplot(1,2,1); pcshow(v_oct,(B_oct)); set(gcf,'color','w'); set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);   xlabel('x (\mu m)');  ylabel('y (\mu m)');  zlabel('z (\mu m)'); title('$z_f = 290 \mu m$'); % caxis([ca1,ca2]); 
    subplot(1,2,2); pcshow(v_oct,(B_octz)); set(gcf,'color','w'); colorbar;  set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);  xlabel('x (\mu m)');  ylabel('y (\mu m)');  zlabel('z (\mu m)');  title('$z_f = 750 \mu m$'); % caxis([ca1,ca2]);

end
