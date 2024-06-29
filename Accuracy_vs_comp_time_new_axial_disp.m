clear all; clc;   % Figure(1) is a plot figure: "Time_Accuracy_new". It is joined with the figure of the script "Comp_Rig_MSR_Tay_Tay2_new.m"

laptop = 'C:\Users\andre\Documents\MATLAB/Matlab_drive_scripts/'; % Desktop personal folder
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
addpath(genpath([path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan']));
addpath(genpath([path_str,'Illumination_functions']));
addpath(genpath([path_str,'saved_heavy_files']));

% Set the Matlab current folder where the folder --> "Simulated_OCT_DW_UL_Bscan" was copied
%cd([path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan'])

tic;

% Folders to add to the path
DW_signal_folder = [path_str,'My_papers_code/Fast_OCT_simulations/saved_data'];
str_folder = [DW_signal_folder,'/' ...%dfgag'];%
    'Acc_vs_Time_865_Nscanx_1_Nscany_1_Nv_64000_Zmin_500_Zmax_1e+03_zf_0_alpha_0_0_NL_1_NA_2.mat'];
%     'Acc_vs_Time_1469_Nscanx_1_Nscany_1_Nv_10240_Zmin_1300_Zmax_1.7e+03_zf_0_alpha_0_0_NL_1_NA_1.mat'];
%     'To_do'];
%     'Acc_vs_Time_1469_Nscanx_1_Nscany_1_Nv_10240_Zmin_1300_Zmax_1.7e+03_zf_0_alpha_0_0_NL_1_NA_1.mat'];

if 0; not(isfolder(remote))&(exist(str_folder) == 2)
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
    refind = 1;
    z_f = 0e-6/refind; % Axial coordinate of the focus of the Gaussian beam: 'z_f > 0'  -->  the mirror is placed at the left of the focal plane of the objective lens
    LX = 80e-6;     % Maximum lateral position of the scatterers in x direction (meters)
    LY = 80e-6;     % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 500e-6;   % minimum of the scatterers axial position (meters)
    maxZ = 1000e-6;  % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 5; % Number of scatterers for ( 10 micron)^3
    Nv = round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    dalpha = 11e-4;
    alpha_min = 0e-3;
    alpha_max = 0e-3;
    NA_tot = 1;
    Nalpha_conf = 2;NA_tot*2 + 1;
    Nalpha_layers = 1; % Number of uniform layers for each combination
    alpha_VEC = 0; linspace(alpha_min + 0.*dalpha/2,alpha_max - 0.*dalpha/2,Nalpha_conf).';% alpha_min + rand(Nalpha_conf,Nalpha_layers) * (alpha_max-alpha_min);
    %   ALPHA = ones(Nalpha_conf, Nalpha_layers).*alpha_VEC;% + (rand(Nalpha_conf, Nalpha_layers)-0.5).*dalpha;
    ALPHA = zeros(Nalpha_conf, Nalpha_layers);
    ALPHA = [0;1e-2];

    %     ALPHA(2:NA_tot+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(2:NA_tot+1,2) = ALPHA(2:NA_tot+1,1);
    %     ALPHA(NA_tot+2:NA_tot*2+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(NA_tot+2:NA_tot*2+1,2) = 1e-3.*10;
    Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
    borders = zeros(1,Nbord+2);        % Borders of the layers of strains
    borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
    borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
    borders(2:Nbord+1)  = minZ + [1000]*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
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
    dscan = 0.7e-6;
    Max_scan = LX - 10e-6; % Maximum lateral scanning in positive and negative x-direction (meters)
    Nscanx = 1; ceil(2*Max_scan/dscan); % Number of A-scan for each B-scan
    if (Nscanx == 1); x_scan_vec = 0e-6; y_scan_vec = -0e-6; Nscany = Nscanx; else; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx); y_scan_vec = 0; Nscany = length(y_scan_vec); end;
    % Nscany = 1; y_scan_vec = 0; Nscanx0 = ceil(Nscanx/2)+1; Nscanx = 2*Nscanx0-1; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx);
    % Uniform refractive index of the medium
    y_scan_vec = 0; x_scan_vec; Nscany = length(y_scan_vec);

    % The following physical parameters have been set as the OCT at UCL
    W_0 = 4.6e-6; % Waist radius of the Gaussian beam at the aperture plane between the collimator and the objective lens
    NA = 0.0972;   % Numerical aperture (Ra = NA*foc2 --> NA = Ra/foc2)
    %   Ra = 3.5e-3; % Radius of the aperture
    %   NA = Ra/foc2;
    lam_min = 1170.5e-9;  % minmum wavelength (meters)
    lam_max = 1407.8e-9;  % maximum wavelength (meters)

    foc1 = 25*1e-3; % focal length of the collimator (meters)
    foc2 = 36*1e-3; % focal length of the objective lens (meters)

    % Dependent parameters
    rho = rand(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)

    %%%%%%%%%%%%%%    Unloaded scatterer region setting    %%%%%%%%%%%%%%%%%%%%%%
    v_U = zeros(Nv,3); % Array of scatterers initialization: 2D array of the scatterers in the unloaded case. This is an array of dimension (Nv,3), where each row consists of 3 numbers that are the (x,y,z) coordinates of each scatterer
    v_U(:,1) = -LX + rand(Nv,1)*2*LX; % (x coordinates of all scatterers)
    v_U(:,2) = -LY + rand(Nv,1)*2*LY; % (y coordinates of all scatterers)
    v_U(:,3) = sort(minZ + rand(Nv,1)*DZ);  % (z coordinates of all scatterers)
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % Display of the chosen parameters (This is useful for simulations with Zeus or Mnemosine)
    if isfolder(remote)
        fprintf(['\n\nParameters:\n Nv = ', num2str(Nv), ' --> Total number of scatterers \n' ...
            ' Nscanx = ', num2str(Nscanx), ' --> Number of A-scans for each B-scan\n Max_scan = ', num2str(Max_scan*1e6) 'mum --> Maximum lateral scanning (left & right)\n LX = ', num2str(LX*1e6), ...
            'mum --> Maximum lateral distance of the scatterers in x direction (meters)\n LY = ', num2str(LY*1e6)...
            ,'mum --> Maximum lateral distance of the scatterers in y direction (meters)\n Nk = ', num2str(Nk),' --> Number of frequencies calculated in the simulation\n rhoMax = ', num2str(rhoMax), ' --> Maximum reflectivity profile of the scatterers\n z_f = ', num2str(z_f*1e6), 'mum --> Axial coordinate of the focus of the Gaussian beam ' ...
            '\n minZ = ', num2str(minZ*1e6), 'mum --> Minimum axial coordinate of the scatterers:\n maxZ = ', num2str(maxZ*1e6), 'mum --> Maximum axial coordinate of the scatterers:\n'...
            , ' Nalpha_layers = ', num2str(Nalpha_layers), '\n Nalpha_conf = ', num2str(Nalpha_conf), ' --> Number of configurations of series of strains. (Each configuration is related to a C-scan)\n alpha_min = ', num2str(alpha_min), '\n alpha_max = ', num2str(alpha_max), '\n refind = ', num2str(refind),'  --> Refractive index (uniform)\n'...
            ,'\nTotal number of A-scans to simulate: ', num2str(Nscanx*Nscany*Nalpha_conf),'\n']);

        answer = input(['\nProceed with those parameters? y/n:\n'],'s');
        while not((answer == "n")|(answer == "y"))
            answer = input(['Answer y or n:\n'],'s');
        end
        if answer == "n"
            error(["Update the parameters list run the script again!"]);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scattered field and modal coefficients calculation (First Born approximation)

    global T_time LB NB Nderx_glob Nderz_glob ak_glob;

    ak_glob = [];
    % subset of data for the test
    iaconf0 = 1; Nk_sub = Nk; min(11,Nk); ik_sub = round(linspace(1,Nk,Nk_sub)).'; % subnumber of wavelenghts
    Nx_sub = Nscanx; min(2,Nscanx); ix_sub = round(linspace(1,Nscanx,Nx_sub)).'; % subnumber of lateral x-scanning
    Ny_sub = Nscany; min(2,Nscany); iy_sub = round(linspace(1,Nscany,Ny_sub)).'; % subnumber of lateral y-scanning

    % Rigorous calculation of the cross-correlation auto-correl. terms (...etc )for  the chosen sample of data.
    T_time = 0; LB = 1; NB = 1; 
    [alpha_output_tmp,IT_CCg_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t(ik_sub),refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    s = seconds(T_time); s.Format = 'hh:mm:ss'; %  T_taylor(iderx,iderz) = T_time;

    z_LP = load_scat_disp(ALPHA(2),BORD(1,:),v_U(:,3)); v_L = v_U; v_L(:,3) = z_LP;
    [alpha_output_tmp,IT_CC_disp] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_L,rho,f_vec_t(ik_sub),refind,W_0/b12,NA,z_f,BORD(1,:),0,x_scan_vec,y_scan_vec);

    T_time = 0; LB = 1; NB = 1;
    [IT_CC_rig, IT_CC_AC_rig,IT_DC_sc_rig,IT_DC_ref_rig,ak_sc_rig,ak_ref_rig] = I_DW_Cscan(v_L,rho,f_vec_t(ik_sub),refind,W_0/b12,NA,z_f,x_scan_vec(ix_sub),y_scan_vec(iy_sub));
    t_rig = T_time; s_rig = seconds(T_time); s_rig.Format = 'hh:mm:ss'; s_rig

    vec_der_x = 4+ [0:4:22]; Nderx = length(vec_der_x);
    vec_der_z = 0:2:10; 0:5; Nderz = length(vec_der_z);
    T_taylor = zeros(Nderx,Nderz); ERROR_T = zeros(Nderx,Nderz);
    IT_CC_TOT = zeros(Nk,Nscanx,Nderx,Nderz); IT_CC_rig = squeeze(IT_CC_rig);

    Nzx = 5; Err_tayg = zeros(Nzx,1); Err_tayg2 = zeros(Nzx,1); Err_tay = zeros(Nzx,1);
    T_tay = zeros(Nzx,1); T_tayg = zeros(Nzx,1); T_tayg2 = zeros(Nzx,1);
    der_xy = zeros(Nzx,2);
    for izx = 1:Nzx
        T_time = 0; LB = 1; NB = 1; Nderx_glob =  vec_der_x(Nderx-Nzx+izx); Nderz_glob = vec_der_z(Nderz-Nzx+izx); der_xy(izx,:) = [Nderx_glob,Nderz_glob];

        [alpha_output,IT_CCg, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
        s = seconds(T_time); s.Format = 'hh:mm:ss'; %  T_taylor(iderx,iderz) = T_time;
        T_tayg(izx) = T_time; 

        fprintf(['Done ', num2str(izx),'/',num2str(Nzx),'\n']);
        Err_tayg(izx) = err_n(IT_CCg(:,1,1,2),IT_CC_rig);
    end
 
    Err_tayg(5)

    % The following lines concerns the file name where I save the data in a struct file.
    str_a = [num2str(alpha_min),'_',num2str(alpha_max),'_NL_',num2str(Nalpha_layers),'_NA_',num2str(Nalpha_conf)];
    str = [DW_signal_folder,'/Acc_vs_Time_',num2str(Nk),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'_alpha_',str_a,'.mat'];
    save(str,'f_vec','f_vec_t','k_vec','Nscanx','Nscany','dscan','Nk','Nv' ,'IT_CC_AC','borders'...
        ,'W_0','ALPHA','NA','Nalpha_conf','b12', 'refind','z_f','LX','LY','minZ','maxZ','v_U','rho',...
        'alpha_output','x_scan_vec','y_scan_vec','IT_DC_sc','IT_DC_ref','akT_sc','akT_ref','Nalpha_layers','alpha_min','alpha_max',...
        'IT_CC_rig', 'T_taylor', 'Err_tay','Err_tayg', 'T_tay','T_tayg2','Err_tayg2',  'T_tayg', 'IT_CC_TOT','ERROR_T', 'Nderx', 'Nderz','IT_CC_AC_rig' ,'IT_DC_sc_rig','IT_DC_ref_rig' );
    disp(str);
 end

if not(isfolder(remote))

    Threshold = 3e-4.*ones(length(T_tayg),1);
    figure(1); clf; NF = 17; NL = 5; Nt1 = 10^(-4.1); Nt2 = 10.^( -0);
    subplot(1,2,1); semilogy(T_tayg,Err_tayg,' k*','LineWidth',NL); hold on; semilogy(T_tayg,Err_tayg,' k--','LineWidth',NL*0.3); ylim([Nt1,Nt2]); xlabel('seconds'); ylabel('Error');   set(gca,'FontSize',NF);  semilogy(T_tayg,Threshold,' r--','LineWidth',NL); title(['Time vs Error']);
end



