% Simulation of a Large OCT C-scan with the Taylor expansion function "15"

clear all; clc;   % Inizialization

TEST = 1;

desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
laptop = 'C:\Users\andre\Documents\MATLAB\Matlab_drive_scripts\'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus

if isfolder(desktop)
    path_str = desktop;
elseif isfolder(remote)
    path_str = remote;
elseif isfolder(laptop)
    path_str = laptop;
end

rmpath(genpath(path_str));
addpath(genpath([path_str,'\My_papers_code\Fast_OCT_simulations\git_hub\']));

% Set the Matlab current folder where the folder --> "Simulated_OCT_DW_UL_Bscan" was copied
% cd([path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan'])

tic;
% Folders to add to the path
saved_data = [path_str,'\My_papers_code\Fast_OCT_simulations\git_hub\saved_data'];


folder_data = 'Several_small_Cscans_Nk_35_Nscanx_29_Nscany_29_Nv_1280_Zmin_0_Zmax_40_zf_0_alpha_0_0.1_NL_1_NA_18.mat';
folder_data = 'Several_small_Cscans_Nk_78_Nscanx_3_Nscany_3_Nv_7200_Zmin_0_Zmax_90_zf_0_alpha_0_0.01_NL_1_NA_3.mat';

str_folder = [saved_data,'/',folder_data];

if   not(isfolder(remote))&(exist(str_folder) == 2)
    DL = load(str_folder); if (isfield(DL,"DL")); DL = DL.DL; end; unpackStruct(DL);
else

    c = 2.997924580105029e+08;  % speed of light

    %%%%%%%%%%%%%%%%%%%     Parameters of the simulation   %%%%%%%%%%%%%%%%%%%%%%%%
    % Scatterers region
    refind = 1;
    z_f = 00e-6/refind; % Axial coordinate of the focus of the Gaussian beam: 'z_f > 0'  -->  the mirror is placed at the left of the focal plane of the objective lens
    LX = 40e-6; % Maximum lateral position of the scatterers in x direction (meters)
    LY = 40e-6; % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 000e-6; % minimum of the scatterers axial position (meters)
    maxZ = 90e-6; % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 6; % Number of scatterers for ( 10 micron)^3
    Nv = 7200; 1270; round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    dalpha = 11e-4;
    alpha_min = 0e-3;
    alpha_max = 1e-2;
    NA_tot = 1;
    Nalpha_conf = 1e3; 1000; NA_tot*2+1; % 
    Nalpha_layers = 1; % Number of uniform layers for each combination
    alpha_VEC =   linspace(alpha_min + 0.*dalpha/2,alpha_max - 0.*dalpha/2,Nalpha_conf).';% alpha_min + rand(Nalpha_conf,Nalpha_layers) * (alpha_max-alpha_min);
    ALPHA = ones(Nalpha_conf, Nalpha_layers).*alpha_VEC;% + (rand(Nalpha_conf, Nalpha_layers)-0.5).*dalpha;
    %ALPHA = zeros(Nalpha_conf, Nalpha_layers);

    %     ALPHA(2:NA_tot+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(2:NA_tot+1,2) = ALPHA(2:NA_tot+1,1);
    %     ALPHA(NA_tot+2:NA_tot*2+1,1) = 1e-3.*([1:NA_tot].*3).'; ALPHA(NA_tot+2:NA_tot*2+1,2) = 1e-3.*10;
    Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
    borders = zeros(1,Nbord+2);        % Borders of the layers of strains
    borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
    borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
    borders(2:Nbord+1)  = minZ +  []*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
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
    dscan = 1.7e-6;
    Max_scan = LX;%-11e-6; LX - 30e-6; % Maximum lateral scanning in positive and negative x-direction (meters)
    Nscanx = ceil(2*Max_scan/dscan); % Number of A-scan for each B-scan  

    if (Nscanx == 1); x_scan_vec = 0e-6; y_scan_vec = -0e-6; Nscany = Nscanx; else;  x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx); y_scan_vec = 0; Nscany = length(y_scan_vec); end;
    % Nscany = 1; y_scan_vec = 0; Nscanx0 = ceil(Nscanx/2)+1; Nscanx = 2*Nscanx0-1; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx);
    % Uniform refractive index of the medium
    y_scan_vec = x_scan_vec; Nscany = length(y_scan_vec);

    % The following physical parameters have been set as the OCT at UCL
    W_0 = 4.6e-6; % Waist radius of the Gaussian beam at the aperture plane between the collimator and the objective lens
    NA = 0.0972;   % Numerical aperture (Ra = NA*foc2 --> NA = Ra/foc2)
    %   Ra = 3.5e-3; % Radius of the aperture
    %   NA = Ra/foc2;
    lam_min = 1170.5e-9;  % minmum wavelength (meters)
    lam_max = 1407.8e-9;  % maximum wavelength (meters)

    foc1 = 25*1e-3; % focal length of the collimator (meters)
    foc2 = 36*1e-3; % focal length of the objective lens (meters)
    %     Ra = 3.5e-3;

    % Dependent parameters
    rho = rand(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%    Unloaded scatterer region setting    %%%%%%%%%%%%%%%%%%%%%%
    v_U = zeros(Nv,3); % Array of scatterers initialization: 2D array of the scatterers in the unloaded case. This is an array of dimension (Nv,3), where each row consists of 3 numbers that are the (x,y,z) coordinates of each scatterer
    v_U(:,1) =   -LX + rand(Nv,1)*2*LX; % (x coordinates of all scatterers)
    v_U(:,2) =    -LY + rand(Nv,1)*2*LY; % (y coordinates of all scatterers)
    v_U(:,3) =   sort(minZ + rand(Nv,1)*DZ);  % (z coordinates of all scatterers)
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

    global T_time LB NB;
    %x_scan_vec2 = linspace(x_scan_vec(1),x_scan_vec(end),301);
    T_time = 0; LB = 1; NB = 1;
    [alpha_output,IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    s = seconds(T_time); s.Format = 'hh:mm:ss'; fprintf(['\nMax error s = ',num2str(max(ERROR_mat)), ', Computational time = ']); s

%     T_time = 0; LB = 1; NB = 1; 
%     [alpha_output2,IT_CC2, IT_CC_AC2,IT_DC_sc2,IT_DC_ref2,akT_sc2,akT_ref2,ERROR_mat2] = Series_Loaded_DW_Cscan_uniform_Taylor16(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD ,ALPHA,x_scan_vec,y_scan_vec);
%     s2 = seconds(T_time); s2.Format = 'hh:mm:ss'; fprintf(['\nMax error s2 = ',num2str(max(ERROR_mat2)), ', Computational time = ']); s2
%     T_time = 0; LB = 1; NB = 1;
%     [IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref] = I_DW_Cscan(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,x_scan_vec,y_scan_vec);

    ik_sub = 1:Nk; ix_sub = ceil(Nscanx*0.5); iy_sub = ceil(Nscany*0.5); 
    T_time = 0; LB = 1; NB = 1; global ak_glob; ak_glob = squeeze(akT_ref(ik_sub,1,1));
    [IT_CC_rig, IT_CC_AC_rig,IT_DC_sc_rig,IT_DC_ref_rig,ak_sc_rig,ak_ref_rig] = I_DW_Cscan(v_U,rho,f_vec_t(ik_sub),refind,W_0/b12,NA,z_f,x_scan_vec(ix_sub),y_scan_vec(iy_sub));
    t_rig = T_time; s_rig = seconds(T_time); s_rig.Format = 'hh:mm:ss'; s_rig
    s_rig_single_A_scan = s_rig;

    z_U = v_U(:,3);
    z_L = load_scat_disp(alpha_max,borders,z_U);
    v_L = v_U; v_L(:,3) = z_L;

%     T_time = 0; LB = 1; NB = 1;
%     [alpha_output2,IT_CC2, IT_CC_AC2,IT_DC_sc2,IT_DC_ref2,akT_sc2,akT_ref2,ERROR_mat2] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_L,rho,f_vec_t,refind,W_0/b12,NA,z_f,borders,ALPHA,x_scan_vec,y_scan_vec);
%     s7 = seconds(T_time); s7.Format = 'hh:mm:ss';
%     %     fprintf(['\nMax error s7 = ',num2str(max(ERROR_mat7)), ', Computational time = ']); s7
% 
%     err_n(akT_sc2,squeeze(akT_sc(:,:,:,end)))

    %     T_time = 0; LB = 1; NB = 1;
    %     [alpha_output7,IT_CC7, IT_CC_AC7,IT_DC_sc7,IT_DC_ref7,akT_sc7,akT_ref7,ERROR_mat7] = Series_Loaded_DW_Cscan_uniform_Taylor5(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    %     s7 = seconds(T_time); s7.Format = 'hh:mm:ss';
    %     fprintf(['\nMax error s7 = ',num2str(max(ERROR_mat7)), ', Computational time = ']); s7
    %     T_time = 0; LB = 1; NB = 1;
    %     [alpha_output8,IT_CC8, IT_CC_AC8,IT_DC_sc8,IT_DC_ref8,akT_sc8,akT_ref8,ERROR_mat8] = Series_Loaded_DW_Cscan_uniform_Taylor11(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    %     s8 = seconds(T_time); s8.Format = 'hh:mm:ss';
    %     fprintf(['\nMax error s8 = ',num2str(max(ERROR_mat8)), ', Computational time = ']); s8

    %     disp(' ');
    %     ERR =  [err_n(IT_CC5,IT_CC7),err_n(IT_CC8,IT_CC7)]

    %     fprintf(['\nComputational time s4 = ']); s4

    %     ERR =  [err_n(IT_CC3,IT_CC),err_n(IT_CC_AC3,IT_CC_AC)]

    % [alpha_output,IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref] = Series_Loaded_DW_Cscan_uniform_Taylor1(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,borders,alpha_vec,x_scan_vec,y_scan_vec);
    % s = seconds(T_time); s.Format = 'hh:mm:ss'; simulation_time = s
    alpha_uniform_vec = ALPHA; %ALPHA = alpha_uniform_vec ;
    simulation_time = s;
    str_a = [num2str(alpha_min),'_',num2str(alpha_max),'_NL_',num2str(Nalpha_layers),'_NA_',num2str(Nalpha_conf)];
    str = [saved_data,'/Several_small_Cscans_Nk_',num2str(Nk),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'_alpha_',str_a,'.mat'];

    wIT = whos('IT_CC'); bIT = wIT.bytes;
    if bIT < 7e8
        save(str,'f_vec','f_vec_t','k_vec','Nscanx','Nscany','dscan','Nk','Nv','IT_CC','IT_CC_AC','borders'...
        ,'W_0','ALPHA','NA','Nalpha_conf','b12','simulation_time','s_rig_single_A_scan','refind','z_f','LX','LY','minZ','maxZ','v_U','rho',...
        'alpha_output','x_scan_vec','y_scan_vec','IT_DC_sc','IT_DC_ref','akT_sc','akT_ref','Nalpha_layers','alpha_min','alpha_max');
    else
         Nalpha_conf_sub = 51; IK_sub = round(linspace(1,Nalpha_conf,Nalpha_conf_sub)); ALPHA_sub = ALPHA(IK_sub,:);
         IT_CC_sub = IT_CC(:,:,:,IK_sub); IT_CC_AC_sub = IT_CC_AC(:,:,:,IK_sub); 
         IT_DC_sc_sub = IT_DC_sc(:,:,:,IK_sub); IT_DC_ref_sub = IT_DC_ref(:,:,:,IK_sub); akT_sc_sub = akT_sc(:,:,:,IK_sub); akT_ref_sub = akT_ref(:,:,:,IK_sub);
         save(str,'Nalpha_conf_sub','ALPHA_sub','f_vec','f_vec_t','k_vec','Nscanx','Nscany','dscan','Nk','Nv','IT_CC_sub','IT_CC_AC_sub','borders'...
        ,'W_0','ALPHA','NA','Nalpha_conf','b12','simulation_time','refind','z_f','LX','LY','minZ','maxZ','v_U','rho',...
        'alpha_output','x_scan_vec','y_scan_vec','IT_DC_sc_sub','s_rig_single_A_scan','IT_DC_ref_sub','akT_sc_sub','akT_ref_sub','Nalpha_layers','alpha_min','alpha_max');
    end
    disp(str)

end

if not(isfolder(remote))
if exist('s_rig_single_A_scan') ;
    Estimated_rigorous_time = Nscanx*Nscany*Nalpha_conf*s_rig_single_A_scan;
    disp(['Estimated time for the rigorous calculation: ']); Estimated_rigorous_time.Format = 'dd:hh:mm:ss'
    disp(['Time for the Taylor calculation: ']); simulation_time
end
    f0 = mean(f_vec); f_max = max(f_vec); f_min = min(f_vec); Df = f_max-f_min;
    cut_off =  1e-1; Wk = 2*pi*(max(f_vec)-f0)/sqrt(log(1/cut_off)); Sk = hanning(Nk); 1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec(:)-f0)./Wk).^2);
    Nyq_z = Nk/Df/4; zmin = 0;-Nyq_z; zmax = maxZ*refind + 50e-6; Nyq_z; Nz = 301; Nx0 = ceil(Nscanx/2);
    [A_tmp2,zbb2] = gen_fft2(-1,f_vec(:)-f0,Sk,2*[zmin,zmax],Nz); Nz0 = min(find(zbb2 >= 0));   zb2 = zbb2 -zbb2(Nz0); zbar = zb2/2;
%     IZ_min = min(find(zbar > 5e-6)); IZ_max = max(find(zbar < maxZ-5e-6)); IZ = IZ_min:IZ_max; Nz2 = length(IZ);
    IZ_min = min(find(zbar > 0e-6)); IZ_max = max(find(zbar < 25e-6)); IZ = IZ_min:IZ_max; Nz2 = length(IZ);
    Nv_tot = Nz2*Nscanx*Nscany;

    if exist("Nalpha_conf_sub")
        IT_CC = IT_CC_sub;
        Nalpha_conf = Nalpha_conf_sub;
        ALPHA = ALPHA_sub;
        IT_CC = IT_CC_sub; 
        clear Nalpha_conf_sub;
    end

    Nalpha_conf_sub_tmp = 9;
    dalpha_sub = ceil(Nalpha_conf./Nalpha_conf_sub_tmp);
    iva_sub = 1:dalpha_sub:Nalpha_conf;
    alpha_conf_sub = ALPHA(iva_sub); Nalpha_conf_sub = length(alpha_conf_sub);

    Nalpha_conf1 = ceil(sqrt(Nalpha_conf_sub));
    Nalpha_conf2 = Nalpha_conf1;

    for ia = 1:Nalpha_conf_sub
        ia
        ia2 = iva_sub(ia);
        BT_CCa = squeeze(gen_fft2(-1,f_vec,IT_CC(:,:,:,ia2).*Sk,zb2));
        v_oct = zeros(Nv_tot,3); B_oct = zeros(Nv_tot,1); MA = max(abs(squeeze(BT_CCa(IZ,:,:))),[],'all');

        for ix = 1:Nscanx
            for iy = 1:Nscany
                it1 = 1 + (iy-1)*Nz2 + (ix-1)*Nz2*Nscanx; it2 = Nz2 + (iy-1)*Nz2 + (ix-1)*Nz2*Nscanx;
                v_oct(it1:it2,:) = [y_scan_vec(ix)*ones(Nz2,1),x_scan_vec(iy)*ones(Nz2,1),zbar(IZ)']*1e6;
                B_oct(it1:it2) = abs(squeeze(BT_CCa(IZ,ix,iy)))./MA;
            end
        end

        eval(['ptCloud',num2str(ia),' = pointCloud(v_oct); ptCloud',num2str(ia),'.Intensity = B_oct;']);
    end

    for ia = 1:Nalpha_conf_sub
        ia
        if Nalpha_conf_sub < 17
            figure(2); hold on; NF = 11;
            subplot(Nalpha_conf1,Nalpha_conf2,ia); eval(['pcshow(ptCloud',num2str(ia),');']);   set(gcf,'color','w'); %set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);
            % title(['\alpha = ', num2str(alpha_conf_sub(ia)*1e3,2), ' m\epsilon']); subplot(1,2,2); pcshow(ptCloud); set(gcf,'color','w'); set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);
        end
    end

end





