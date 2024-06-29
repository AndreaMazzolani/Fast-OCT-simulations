% Figure(1) is a plot of the figure: "Time_Accuracy_new". It is joined from the figure of the script "Accuracy_vs_comp_time_new.m"

clear all; clc;   % Inizialization
rng(1);
TEST = 1;

desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus
laptop = 'C:/Users/andre/Documents/MATLAB/Matlab_drive_scripts/'; % personal folder in my Laptop (October 2022)

if isfolder(desktop)
    path_str = desktop;
elseif isfolder(laptop)
    path_str = laptop;
elseif isfolder(remote)
    path_str = remote;
else error(['The path-folder does not exist!']);
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
    'Test_dataset_new_Nkm2_100_Nkm_2048_NkTmax_200_Nvm_200000_Nvm2_10000_NmT_200_Nind_k_10_large.mat'];
%    'Test_dataset_new_Nkm2_100_Nkm_2048_NkTmax_200_Nvm_200000.mat']; % Good  from Clio
%  'Test_dataset_new_Nkm2_100_Nkm_2048_NkTmax_200_Nvm_200000_Nvm2_10000_NmT_200_Nind_k_10.mat']; % Good from Mnem
%     'Test_dataset_new_Nkm2_100_Nkm_2048_NkTmax_200_Nvm_200000_Nvm2_10000_NmT_200_Nind_k_10.mat'];
%     'Test_dataset_new_Nkm2_1000_Nkm_148_NkTmax_5_Nvm_2000_Nvm2_10000_NmT_5_Nind_k_10.mat'];
%      'Test_dataset_new_'];

if  not(isfolder(remote))&(exist(str_folder) == 2)
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
    LX = 100e-6;     % Maximum lateral position of the scatterers in x direction (meters)
    LY = 100e-6;     % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 1200e-6;   % minimum of the scatterers axial position (meters)
    maxZ =  2210e-6;  % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 5; % Number of scatterers for ( 10 micron)^3
    Nv = 10000; round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    dalpha = 11e-4;
    alpha_min = -5e-3;
    alpha_max = 5e-3;
    NA_tot = 1;
    Nalpha_conf = 11;NA_tot*2+1;
    Nalpha_conf2 = 31;NA_tot*2+1; Nalpha_conf = Nalpha_conf2;
    Nalpha_layers = 1; % Number of uniform layers for each combination
    alpha_VEC = linspace(alpha_min + 0.*dalpha/2,alpha_max - 0.*dalpha/2,Nalpha_conf).';% alpha_min + rand(Nalpha_conf,Nalpha_layers) * (alpha_max-alpha_min);
    ALPHA = ones(Nalpha_conf, Nalpha_layers).*alpha_VEC;% + (rand(Nalpha_conf, Nalpha_layers)-0.5).*dalpha;
%     ALPHA = zeros(Nalpha_conf, Nalpha_layers);
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

    % Lateral scanning p-arameters
    dscan = 0.7e-6; Max_scan = LX - 10e-6; % Maximum lateral scanning in positive and negative x-direction (meters)
    Nscanx = 1; ceil(2*Max_scan/dscan); % Number of A-scan for each B-scan
    if (Nscanx == 1); x_scan_vec = 0e-6; y_scan_vec = -0e-6; Nscany = Nscanx; else; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx); y_scan_vec = 0; Nscany = length(y_scan_vec); end;
    y_scan_vec = 0; x_scan_vec; Nscany = length(y_scan_vec);

    % The following physical parameters have been set as the OCT at UCL
    W_0 = 4.6e-6; % Waist radius of the Gaussian beam at the aperture plane between the collimator and the objective lens
    NA = 0.0972;   % Numerical aperture (Ra = NA*foc2 --> NA = Ra/foc2)
    foc1 = 25*1e-3; foc2 = 36*1e-3; % focal length of the objective lens (meters)
    b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)

    %%%%%%%%%%%%%%    Unloaded scatterer region setting    %%%%%%%%%%%%%%%%%%%%%%
    v_U = zeros(Nv,3); % Array of scatterers initialization: 2D array of the scatterers in the unloaded case. This is an array of dimension (Nv,3), where each row consists of 3 numbers that are the (x,y,z) coordinates of each scatterer
    v_U(:,1) = -LX  + rand(Nv,1)*2*LX; % (x coordinates of all scatterers)
    v_U(:,2) = -LY + rand(Nv,1)*2*LY; % (y coordinates of all scatterers)
    v_U(:,3) = sort(minZ + rand(Nv,1)*DZ);  % (z coordinates of all scatterers)
    rhoMax = 3e-4; rho = rand(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%     Spatial frequency setting  'f = 1/lambda'    %%%%%%%
    lam_min = 1170.5e-9; lam_max = 1407.8e-9;  % maximum wavelength (meters)
    f_min = 1./lam_max; f_max = 1./lam_min;  % maximum frequency of the simulated spectrum
    f0 = ( f_max + f_min)/2; Df = f_max - f_min; % frequencies spectrum width
    Nk_nyq = 4*max(abs(v_U(:,3).*refind))*Df; % Minimum number of frequencies able to simulate the deepest axial position (coming from the Nyquist frequency and rates relation)
    Nk = ceil(Nk_nyq*1.5); % Number of wavenumbers, they are set as 1.5 times 'Nk_nyq' ( A safe number greater than 'Nk_nyq')
    f_vec = linspace(f_min,f_max,Nk); f_vec_t = f_vec.*c; k_vec = 2*pi*f_vec; % Arrays of wavenumbers;

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
  
    rho_max = max(sqrt(v_U(:,1).^2 + v_U(:,2).^2));
    D_vec = linspace(0,100e-6,1001); Upper_bound_B = rho_max*NA + maxZ*NA^2/refind + refind*D_vec/2;
    % clf; plot(D_vec*1e6,Upper_bound_B)

    Nkm2 = 50; Nkm = 2048;  NkTmax = 200;  Nvm = 100000; Nvm2 = 10000; NmT = 200;
%     Nkm2 = 1000; Nkm = 148; Nvm = 2000; Nvm2 = 10000; NmT = 5; NkTmax = 5; 
    
    Nind_k = 10; ind_vec_k = ceil([1:Nind_k]*Nkm/Nind_k); 
    str_save = [str_folder,'Nkm2_',num2str(Nkm2),'_Nkm_',num2str(Nkm),'_NkTmax_',num2str(NkTmax),...
        '_Nvm_',num2str(Nvm),'_Nvm2_',num2str(Nvm2),'_NmT_',num2str(NmT),'_Nind_k_',num2str(Nind_k),'.mat']
    list_save = {'Nind_k';'Nkm2';'Nkm';'Nvm';'Nvm2';'NmT';'NkTmax';'f_vec';'Nscanx';'Nscany';'Nk';'minZ';'maxZ';'x_scan_vec';'y_scan_vec';'alpha_min';'alpha_max'}; Nlist = length(list_save);
    for il = 1:Nlist; eval(['DL.', list_save{il},' = ',list_save{il},';']); end;  save(str_save,'DL');

    Nk_vect = ceil(Nkm2.*[1:NkTmax]); Nv_vect = ceil(Nvm2.*[1:NmT]); 
    T_taylor_Nv = zeros(NmT,1); T_taylor_Nv_rig = zeros(NmT,1); T_taylor_Nv_no_T2 = zeros(NmT,Nalpha_conf); T_taylor_Nv_no_T_T2 = zeros(NmT,Nalpha_conf); T_taylor_Nv_no_MSR_no_T2 = zeros(NmT,Nalpha_conf);
    err_Nv = zeros(NmT,1); err_Nv_no_T2 = zeros(NmT,1); err_Nv_no_T_T2 = zeros(NmT,1); err_Nv_no_MSR_no_T2 = zeros(NmT,1);
    T_taylor_Nk = zeros(NkTmax,1); T_taylor_Nk_rig = zeros(NkTmax,1); T_taylor_Nk_no_T2 = zeros(NkTmax,Nalpha_conf); T_taylor_Nk_no_T_T2 = zeros(NkTmax,Nalpha_conf); T_taylor_Nk_no_MSR_no_T2 = zeros(NkTmax,Nalpha_conf);
    err_Nk = zeros(NkTmax,1); err_Nk_no_T2 = zeros(NkTmax,1); err_Nk_no_T_T2 = zeros(NkTmax,1); err_Nk_no_MSR_no_T2 = zeros(NkTmax,1);

    for im = 1:NmT
        Nv_tot = Nv_vect(im);  
        v_Ut = zeros(Nv_tot,3); v_Ut(:,1) = -LX + rand(Nv_tot,1)*2*LX; v_Ut(:,2) = -LY + rand(Nv_tot,1)*2*LY; % (y coordinates of all scatterers)
        v_Ut(:,3) = sort(minZ + rand(Nv_tot,1)*DZ); rhot = rand(Nv_tot,1).*rhoMax;

        % ak_ref definition
        global T_time LB NB ak_glob; ak_glob = []; T_time = 0; LB = 1; NB = 1; f_vec_t2 = linspace(f_vec_t(1),f_vec_t(Nk),Nkm);
        [alpha_output,IT_CC_test, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Ut(1,:),0,f_vec_t2,refind,W_0/b12,NA,z_f,BORD(1,:),0*ALPHA(1,:),x_scan_vec,y_scan_vec);
        ak_glob = ak_ref;
    
        % Taylor "MSR+T1+T2" (total)
        T_time = 0; LB = 1; NB = 1;% ak_glob = ak_glob(2);
        [alpha_output0,IT_CC] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Ut,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
        T_taylor_Nv(im) = T_time; disp(['Done ',num2str(im),'/',num2str(NmT)]);
%     end
        IT_CC_noT2 = IT_CC.*0; IT_CC_noMSR = IT_CC.*0; IT_CC_noT = IT_CC.*0; 
        
        ia = Nalpha_conf;
%         fprintf('Nv est.: (%d,%d)/(%d,%d)\n', im,ia,NmT,Nalpha_conf);
        % Taylor "MSR+T1" (No T2)
        alpha_vec = ALPHA(ia,:); bord_vec = BORD(ia,:); zL = load_scat_disp(alpha_vec,bord_vec,v_Ut(:,3)); v_Lt = v_Ut; v_Lt(:,3) = zL; T_time = 0; LB = 1; NB = 1; 
        [alpha_output0,IT_CC_noT2_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nv_no_T2(im,ia) = T_time; IT_CC_noT2(:,:,:,ia) = IT_CC_noT2_tmp;
        
        T_time = 0; LB = 1; NB = 1;  % Taylor  "T1" (No MSR-No T2)
        [alpha_output0,IT_CC_noMSR_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor_noMSR(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nv_no_MSR_no_T2(im,ia) = T_time; IT_CC_noMSR(:,:,:,ia) = IT_CC_noMSR_tmp;
        
        T_time = 0; LB = 1; NB = 1;  % Taylor "MSR" (No T1-T2)
        [alpha_output0,IT_CC_noT_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor_noT(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nv_no_T_T2(im,ia) = T_time; IT_CC_noT(:,:,:,ia) = IT_CC_noT_tmp;
        
        T_time = 0; LB = 1; NB = 1; ak_glob = ak_glob(ind_vec_k);
        [IT_CC_rig, IT_CC_AC_rig,IT_DC_sc_rig,IT_DC_ref_rig,ak_sc_rig,ak_ref_rig] = I_DW_Cscan(v_Lt,rhot,f_vec_t2(ind_vec_k),refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);
        T_taylor_Nv_rig(im) = T_time; 
        err_Nv(im) = err_n(IT_CC_rig,IT_CC(ind_vec_k,Nscanx,Nscany,Nalpha_conf));
        err_Nv_no_T2(im) = err_n(IT_CC(:,:,:,ia),IT_CC_noT2(:,:,:,ia));
        err_Nv_no_T2(im) = err_n(IT_CC(:,:,:,ia),IT_CC_noMSR(:,:,:,ia));
        err_Nv_no_T2(im) = err_n(IT_CC(:,:,:,ia),IT_CC_noT(:,:,:,ia));
    end
    list_save2 = {'Nk_vect','Nv_vect','T_taylor_Nv','T_taylor_Nv_rig','T_taylor_Nv_no_T2','T_taylor_Nv_no_T_T2',...
        'T_taylor_Nv_no_MSR_no_T2','err_Nv',  'err_Nv_no_T2', 'err_Nv_no_T_T2' ,'err_Nv_no_MSR_no_T2'}; Nlist2 = length(list_save2);
    for il = 1:Nlist2; eval(['DL.', list_save2{il},' = ',list_save2{il},';']); end;   save(str_save,'DL');
  
    for ik = 1:NkTmax-1;
        v_Ut = zeros(Nvm,3); v_Ut(:,1) = -LX + rand(Nvm,1)*2*LX; v_Ut(:,2) = -LY + rand(Nvm,1)*2*LY; % (y coordinates of all scatterers)
        v_Ut(:,3) = sort(minZ + rand(Nvm,1)*DZ); rhot = rand(Nvm,1).*rhoMax;
        Nk2 = Nk_vect(ik);  f_vec_t2 = linspace(f_vec_t(1),f_vec_t(Nk),Nk2);

        % ak_ref definition
        ak_glob = []; T_time = 0; LB = 1; NB = 1;  
        [alpha_output,IT_CC_test, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Ut(1,:),0,f_vec_t2,refind,W_0/b12,NA,z_f,BORD(1,:),0*ALPHA(1,:),x_scan_vec,y_scan_vec);
        ak_glob = ak_ref;
    
        % Taylor "MSR+T1+T2" (total)
        T_time = 0; LB = 1; NB = 1;% ak_glob = ak_glob(2);
        [alpha_output0,IT_CC] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Ut,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
        T_taylor_Nk(ik) = T_time; disp(['Done ',num2str(ik),'/',num2str(NkTmax)]);
%     end
        IT_CC_noT2 = IT_CC.*0; IT_CC_noMSR = IT_CC.*0; IT_CC_noT = IT_CC.*0; 
        ia = Nalpha_conf;
%         fprintf('Nk est.: (%d,%d)/(%d,%d)\n', ik,ia,NkTmax,Nalpha_conf);
        % Taylor "MSR+T1" (No T2)
        alpha_vec = ALPHA(ia,:); bord_vec = BORD(ia,:); zL = load_scat_disp(alpha_vec,bord_vec,v_Ut(:,3)); v_Lt = v_Ut; v_Lt(:,3) = zL; T_time = 0; LB = 1; NB = 1; 
        [alpha_output0,IT_CC_noT2_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nk_no_T2(ik,ia) = T_time; IT_CC_noT2(:,:,:,ia) = IT_CC_noT2_tmp;
        
        T_time = 0; LB = 1; NB = 1;  % Taylor  "T1" (No MSR-No T2)
        [alpha_output0,IT_CC_noMSR_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor_noMSR(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nk_no_MSR_no_T2(ik,ia) = T_time; IT_CC_noMSR(:,:,:,ia) = IT_CC_noMSR_tmp;
        
        T_time = 0; LB = 1; NB = 1;  % Taylor "MSR" (No T1-T2)
        [alpha_output0,IT_CC_noT_tmp] = Series_Loaded_DW_Cscan_uniform_Taylor_noT(v_Lt,rhot,f_vec_t2,refind,W_0/b12,NA,z_f,[min(zL),max(zL)],0,x_scan_vec,y_scan_vec);
        T_taylor_Nk_no_T_T2(ik,ia) = T_time; IT_CC_noT(:,:,:,ia) = IT_CC_noT_tmp;
    
        ind_vec_k_tmp = [1,Nk2]; Nind_k2 = length(ind_vec_k_tmp); ak_glob = ak_glob(ind_vec_k_tmp);
        if ik == 1; 
            T_time = 0; LB = 1; NB = 1; 
            [IT_CC_rig_tmp, IT_CC_AC_rig,IT_DC_sc_rig,IT_DC_ref_rig,ak_sc_rig,ak_ref_rig] = I_DW_Cscan(v_Lt,rhot,f_vec_t2(ind_vec_k_tmp),refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);
            T_taylor_Nk_rig(ik,ia) = T_time; IT_CC_rig = IT_CC_rig_tmp;
        end
        err_Nk(ik) = err_n(IT_CC_rig,IT_CC(ind_vec_k_tmp,Nscanx,Nscany,Nalpha_conf));
        err_Nk_no_T2(ik) = err_n(IT_CC(:,:,:,ia),IT_CC_noT2(:,:,:,ia));
        err_Nk_no_T2(ik) = err_n(IT_CC(:,:,:,ia),IT_CC_noMSR(:,:,:,ia));
        err_Nk_no_T2(ik) = err_n(IT_CC(:,:,:,ia),IT_CC_noT(:,:,:,ia));
    end
    list_save3 = {'Nind_k2','T_taylor_Nk','T_taylor_Nk_rig','T_taylor_Nk_no_T2','T_taylor_Nk_no_T_T2',...
        'T_taylor_Nk_no_MSR_no_T2','err_Nk',  'err_Nk_no_T2', 'err_Nk_no_T_T2' ,'err_Nk_no_MSR_no_T2'}; Nlist3 = length(list_save3);
    for il = 1:Nlist3; eval(['DL.', list_save3{il},' = ',list_save3{il},';']); end;  
    
%     save(str_save,'DL'); disp(str_save)
    str_save2 = [str_folder,'Nkm2_',num2str(Nkm2),'_Nkm_',num2str(Nkm),'_NkTmax_',num2str(NkTmax),...
        '_Nvm_',num2str(Nvm),'_Nvm2_',num2str(Nvm2),'_NmT_',num2str(NmT),'_Nind_k_',num2str(Nind_k),'_large.mat']
    
    Nalpha_conf2 = 51; DL.Nalpha_conf2 = Nalpha_conf2; save(str_save2,'DL'); disp(str_save2)

    if strcmp(path_str,remote)
        setenv('TOUPLOAD', str_save);
        ! scp -r $TOUPLOAD "amazzolani@128.40.68.241:/home/amazzolani/Matlab_drive_scripts/My_papers_code/Fast_OCT_simulations/saved_data/${TOUPLOAD}"
       % ! scp -r $TOUPLOAD "andrea:/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/My_papers_code/Fast_OCT_simulations/saved_data/"
%        ! scp -r $TOUPLOAD "C:/Users/andre/Documents/MATLAB/Matlab_drive_scripts/My_papers_code/Fast_OCT_simulations/saved_data/${TOUPLOAD}"
    end
    % /usr/local/MATLAB/R2022a/MATLAB_26_04_2022/Matlab_drive_scripts/Simulated_OCT_DW_seriesL_uniform_signals/Unloaded_series_of_Loaded_Bscans_uniform_Debye.mlx amazzolani@128.40.68.241:/home/amazzolani/Matlab_drive_scripts/Simulated_OCT_DW_UL_Bscan
end

if not(isfolder(remote))
    Nalpha_conf = size(T_taylor_Nv_no_T2,2); N_mean = 5;

    if (exist("Nalpha_conf2")); ia = Nalpha_conf; Nalpha_conf = Nalpha_conf2; end;

    secs = 1;
    mins = secs*60; 
    hours = mins*60;
    days = hours*24;
    weeks = days*7;
    months = days*30;
    years = months*12;

    vec_time = [secs mins hours days weeks months];
    vec_str_time = {'1 sec','1 min','1 hour','1 day','1 week','1 month'};
    
    [Tv_rig,ind_v] = slit_mean(T_taylor_Nv_rig./Nind_k*Nkm*Nalpha_conf,N_mean); Tv_MTP = slit_mean(T_taylor_Nv,N_mean); Tv_MT_ = slit_mean(T_taylor_Nv_no_T2(:,ia),N_mean)*Nalpha_conf;  Tv__T_ = slit_mean(T_taylor_Nv_no_MSR_no_T2(:,ia),N_mean)*Nalpha_conf; Tv_M__ = slit_mean(T_taylor_Nv_no_T_T2(:,ia),N_mean)*Nalpha_conf;
    [Tk_rig,ind_k] = slit_mean(T_taylor_Nk_rig(1,ia)./Nind_k2*Nk_vect(:)*Nalpha_conf,N_mean); Tk_MTP = slit_mean(T_taylor_Nk,N_mean); Tk_MT_ = slit_mean(T_taylor_Nk_no_T2(:,ia),N_mean)*Nalpha_conf;  Tk__T_ = slit_mean(T_taylor_Nk_no_MSR_no_T2(:,ia),N_mean)*Nalpha_conf; Tk_M__ = slit_mean(T_taylor_Nk_no_T_T2(:,ia),N_mean)*Nalpha_conf;

    figure(1); clf; NF = 10; NL = 1.8; sgtitle('Computational time of one A-scan with 50 loadings');
    subplot(1,2,1); semilogy(Nv_vect(ind_v),Tv_rig,'k--','Linewidth',NL,'Linewidth',NL); hold on; ylim([0,weeks*11]);    yticks(vec_time); yticklabels(vec_str_time);
    semilogy(Nv_vect(ind_v),Tv__T_,'Linewidth',NL); semilogy(Nv_vect(ind_v),Tv_M__,'Linewidth',NL);  semilogy(Nv_vect(ind_v),Tv_MT_,'Linewidth',NL); semilogy(Nv_vect(ind_v),Tv_MTP,'Linewidth',NL);  title(['Nk = ',num2str(Nkm)],'Linewidth',NL); xlabel('Nv (Num. scatterers)'); set(gca, 'FontSize',NF);  
    subplot(1,2,2); title('Nv = '); semilogy(Nk_vect(ind_k),Tk_rig,'k--','Linewidth',NL); hold on; ylim([0,weeks*11]);     % yticks(vec_time); yticklabels(vec_str_time);
    semilogy(Nk_vect(ind_k),Tk__T_,'Linewidth',NL); semilogy(Nk_vect(ind_k),Tk_M__,'Linewidth',NL); semilogy(Nk_vect(ind_k),Tk_MT_,'Linewidth',NL); semilogy(Nk_vect(ind_k),Tk_MTP,'Linewidth',NL);  legend('rig','T','M','M+T','M+T+A','Linewidth',NL); xlabel('Nk (Num. wavenumbers)'); title(['Nv = ',num2str(Nvm)],'Linewidth',NL); set(gca, 'FontSize',NF,'YTick', []); 
    
end
