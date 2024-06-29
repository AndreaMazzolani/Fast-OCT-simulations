clear all; clc;

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
if isfolder(desktop)
    path_str = desktop;
elseif isfolder(remote)
    path_str = remote;
end

addpath(genpath([path_str,'handle_functions']));
addpath(genpath([path_str,'OCT_models/Simulated_OCT_DW_UL_Bscan']));
addpath(genpath([path_str,'OCE_models/Phase_difference_based_methods']));
addpath(genpath([path_str,'Illumination_functions']));
addpath(genpath([path_str,'hyperion']));
addpath(genpath([path_str,'saved_heavy_files']));
% cd(path_str)

% addpath(genpath([path_str,'Debye_fields_OCT_UCL_Nk_3/']));
% addpath('FDTD_simulations');

DW_signal_folder = [path_str,'saved_heavy_files'];%My_papers_code/Fast_OCT_simulations/saved_data'];
str_folder = [DW_signal_folder,'/' ...
  'Several_small_Cscans_Nk_130_Nscanx_36_Nscany_36_Nv_7500_Zmin_0_Zmax_150_zf_0_alpha_0_0.03_NL_1_NA_1000.mat'];
  %'OCE_comparisonanx_1_Nscany_1_Nv_5000_Zmin_0_Zmax_500_zf_0_alpha_0_0.002_NL_3_NA_2.mat'];
%     'OCE_comparison_432_Nscanx_21_Nscany_1_Nv_9600_Zmin_0_Zmax_500_zf_0_alpha_0_0.002_NL_3_NA_2.mat'];
% 'OCE_comparison_433_Nscanx_31_Nscany_1_Nv_5000_Zmin_0_Zmax_500_zf_0_alpha_0_0.002_NL_3_NA_2.mat'];
c = 2.997924580105029e+08;  % speed of light


if not(isfolder(remote))&(exist(str_folder) == 2)
    DL = load(str_folder); if (isfield(DL,"DL")); DL = DL.DL; end
    unpackStruct(DL);
else

    %%%%%%%%%%%%%%%%%%%     Parameters of the simulation   %%%%%%%%%%%%%%%%%%%%%%%%
    % Scatterers region
    refind = 1;
    z_f = 500e-6/refind; % Axial coordinate of the focus of the Gaussian beam: 'z_f > 0'  -->  the mirror is placed at the left of the focal plane of the objective lens
    LX = 96e-6;     % Maximum lateral position of the scatterers in x direction (meters)
    LY = 96e-6;     % Maximum lateral position of the scatterers in y direction (meters)
    minZ = 000e-6;   % minimum of the scatterers axial position (meters)
    maxZ = 1024e-6;  % maximum of the scatterer axial  position (meters) (meters)
    DZ = maxZ - minZ; % Width of the axial interval
    nn3_dens = 3; % Number of scatterers for ( 10 micron)^3
    Nv = 5e5; round( 2*LX*2*LY*DZ/(10e-6^3)*nn3_dens); % The number of scatterers is usually set as a function of the region volume

    % Layers of strains
    Nalpha_conf = 101;
    Nalpha_layers = 1; % Number of uniform layers for each combination
    %     ALPHA(1,:) = [0];
    %     ALPHA(2,:) = [3]*1e-3;
    alpha_min = -5e-3; alpha_max = 5e-3; 
    ALPHA =  linspace(alpha_min,alpha_max,Nalpha_conf).';
    %   alpha_min = min(ALPHA(:)); alpha_max = max(ALPHA(:));
    Nbord = Nalpha_layers-1; % it defines Nbord + 1 uniform strains for each C-scan
    borders = zeros(1,Nbord+2);        % Borders of the layers of strains
    borders(1) = minZ;                  % 1-st boarder must be smaller or equal to the shallowest scatterer z-coordinate
    borders(Nbord+2) = maxZ;           % the last boarder must be greater or equal to the deepest scatterer z-coordinate
    borders(2:Nbord+1)  = [ ]*1e-6;  % Choose the coordinates of the borders of the middle layers! ( borders (2:N-1));
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
    dscan = 1.5e-6;
    Max_scan = LX - 0e-6; % Maximum lateral scanning in positive and negative x-direction (meters)
    Nscanx = 128; ceil(2*Max_scan/dscan); % Number of A-scan for each B-scan
    if (Nscanx == 1); x_scan_vec = 0e-6; y_scan_vec = -0e-6; Nscany = Nscanx; else; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx); y_scan_vec = 0; Nscany = length(y_scan_vec); end;
    % Nscany = 1; y_scan_vec = 0; Nscanx0 = ceil(Nscanx/2)+1; Nscanx = 2*Nscanx0-1; x_scan_vec = linspace(-Max_scan,Max_scan,Nscanx);
    % Uniform refractive index of the medium
    y_scan_vec =  x_scan_vec; Nscany = length(y_scan_vec);

    % The following physical parameters have been set as the OCT at UCL
    W_0 = 4.6e-6; % Waist radius of the Gaussian beam at the aperture plane between the collimator and the objective lens
    NA = 0.0972;   % Numerical aperture (Ra = NA*foc2 --> NA = Ra/foc2)
    lam_min = 1170.5e-9;  % minmum wavelength (meters)
    lam_max = 1407.8e-9;  % maximum wavelength (meters)

    foc1 = 25*1e-3; foc2 = 36*1e-3;

    % Dependent parameters
    rho = rand(1,Nv).*rhoMax; % Array of the reflectivity profiles of all scatterers (It is set random)
    b12 = foc1/foc2;  % ratio between the focal lenghts ( This parameter is used in the Debye-Wolf integral)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


    global T_time LB NB;
    T_time = 0; LB = 1; NB = 1;
    [alpha_output,IT_CC, IT_CC_AC,IT_DC_sc,IT_DC_ref,akT_sc,akT_ref,ERROR_mat,ak_ref] = Series_Loaded_DW_Cscan_uniform_Taylor15(v_U,rho,f_vec_t,refind,W_0/b12,NA,z_f,BORD,ALPHA,x_scan_vec,y_scan_vec);
    s = seconds(T_time); s.Format = 'hh:mm:ss';   T_taylor = T_time;
    fprintf(['\nMax error s = ',num2str(max(ERROR_mat)), ', Computational time = ']); s


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    OCT and OCE  SIGNAL RECONSTRUCTION FROM RAW DATA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    str_a = [num2str(alpha_min),'_',num2str(alpha_max),'_NL_',num2str(Nalpha_layers),'_NA_',num2str(Nalpha_conf)];

    str1 = [DW_signal_folder,'/OCE_comparison_',num2str(Nk),'_Nscanx_',num2str(Nscanx),'_Nscany_',num2str(Nscany),'_Nv_',num2str(Nv),'_Zmin_',num2str(minZ*1e6),'_Zmax_',num2str(maxZ*1e6,3),'_zf_',num2str(z_f*1e6,3),'_alpha_',str_a,'.mat'];
    save(str1,'f_vec','f_vec_t','k_vec','Nscanx','Nscany','dscan','Nk','Nv','IT_CC','borders','Nbord','BORD'...
        ,'W_0','ALPHA','NA','Nalpha_conf','b12','refind','z_f','LX','LY','minZ','maxZ','v_U','rho',...
        'x_scan_vec','y_scan_vec','Nalpha_layers','alpha_min','alpha_max',   'T_taylor');
    disp(str1)
end

if size(IT_CC,4) > 2; error('For this example only two compressed version of the same OCT signal!'); end

IU_CC = squeeze(IT_CC(:,:,:,1)); IL_CC = squeeze(IT_CC(:,:,:,2));

f0 = mean(f_vec); f_max = max(f_vec);
cut_off = 1e-3; Wk = 2*pi*(f_max-f0)/sqrt(log(1/cut_off)); Sk =  1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec(:)-f0)./Wk).^2);
Df = f_vec(Nk)-f_vec(1); Nyq_z = Nk/Df/4; zmin = 50e-6;-Nyq_z; zmax = 450e-6;Nyq_z; Nz = 1501;
[A_tmp,zb2] = gen_fft2(-1,f_vec(:),IU_CC(:,1),2*[zmin,zmax],Nz);  zbar = zb2/2; %zb2 = zbb2 - zbb2(Nz0); zb2 = zb2(:); zbar = zb2./2;
BU_CC = gen_fft2(-1,f_vec,IU_CC.*Sk,zb2); BL_CC = gen_fft2(-1,f_vec,IL_CC.*Sk,zb2);

% eps1 = 1e-10;
% DBU_CC2 = gen_fft2(-1,f_vec,IU_CC.*Sk,2.*(zbar + eps1)); DBL_CC2 = gen_fft2(-1,f_vec,IL_CC.*Sk,2.*(zbar + eps1));
% DBU_CC1 = gen_fft2(-1,f_vec,IU_CC.*Sk,2.*(zbar - eps1)); DBL_CC1 = gen_fft2(-1,f_vec,IL_CC.*Sk,2.*(zbar - eps1));
% DBU_CC = (DBU_CC2 - DBU_CC1)./2./eps1; DBL_CC = (DBL_CC2 - DBL_CC1)./2./eps1;
%
% DBU_CCr2 = gen_fft2(-1,f_vec,IU_CC_rig.*Sk,2.*(zbar + eps1)); DBL_CCr2 = gen_fft2(-1,f_vec,IL_CC_rig.*Sk,2.*(zbar + eps1));
% DBU_CCr1 = gen_fft2(-1,f_vec,IU_CC_rig.*Sk,2.*(zbar - eps1)); DBL_CCr1 = gen_fft2(-1,f_vec,IL_CC_rig.*Sk,2.*(zbar - eps1));
% DBU_CC_rig = (DBU_CCr2 - DBU_CCr1)./2./eps1; DBL_CC_rig = (DBL_CCr2 - DBL_CCr1)./2./eps1;
%

if not(strcmp(path_str,remote))

    dphi_wrap =  (angle(BU_CC./BL_CC));
    dphi = unwrap(angle(BU_CC./BL_CC));

    % clf; imagesc(dphi_wrap); colorbar; axis equal;

    if length(size(BU_CC)) > 2
        IZ_min = min(find(zbar > 0e-6)); IZ_max = max(find(zbar < 100e-6)); IZ = IZ_min:IZ_max; Nz2 = length(IZ);
        Nv_tot = Nz2*Nscanx*Nscany; v_oct = zeros(Nv_tot,3); B_oct = zeros(Nv_tot,1);
        MA = max( (abs(BU_CC(:))));

        for ix = 1:Nscanx
            ix
            for iy = 1:Nscany
                it1 = 1 + (iy-1)*Nz2 + (ix-1)*Nz2*Nscanx; it2 = Nz2 + (iy-1)*Nz2 + (ix-1)*Nz2*Nscanx;
                v_oct(it1:it2,:) = [x_scan_vec(ix)*ones(Nz2,1),y_scan_vec(iy)*ones(Nz2,1),zbar(IZ)']*1e6;
                B_oct(it1:it2) =  squeeze(dphi_wrap(IZ,ix,iy)); %(abs(BU_CC(IZ,ix,iy)))./MA;
            end
        end

        ptCloud = pointCloud(v_oct); ptCloud.Intensity = B_oct;
        figure(13); clf; NF = 11;
        pcshow(ptCloud); set(gcf,'color','w'); set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]); colorbar; %colormap('gray');

        fit_range = 80e-6; k0 = f0*2*pi; lamb0 = 1./f0;
        dz = zbar(2) - zbar(1); dzphi = dphi / (-2*k0);



        clear strain;
        for iy = 1:Nscany
            iy
            dzphiy = squeeze(dzphi(:,:,iy,1));
            [strainy, zp] = conventional_reconstruct(zbar, dzphiy.', fit_range, 0, true);
            strain(:,iy,:) = strainy;
        end

        IZp_min = min(find(zp > 0e-6)); IZp_max = max(find(zp < 100e-6)); IZp = IZp_min:IZp_max; Nzp2 = length(IZp);
        Nvp_tot = Nzp2*Nscanx*Nscany; %vp_oct = zeros(Nvp_tot,3); Bp_oct = zeros(Nvp_tot,1);

        clear vp_oct Bp_oct;
        for ix = 1:Nscanx
            ix
            if ix > ceil(Nscanx/2)
                IY = 1:Nscany;
            else
                IY =1:Nscany; ceil(Nscany/2):Nscany;
            end

            for iy = IY;%1:Nscany
                it1 = 1 + (iy-1)*Nzp2 + (ix-1)*Nzp2*Nscanx; it2 = Nzp2 + (iy-1)*Nzp2 + (ix-1)*Nzp2*Nscanx;
                vp_oct(it1:it2,:) = [x_scan_vec(ix)*ones(Nzp2,1),y_scan_vec(iy)*ones(Nzp2,1),zp(IZp)']*1e6;
                Bp_oct(it1:it2) =  strain(ix,iy,IZp); %(abs(BU_CC(IZ,ix,iy)))./MA;
            end
        end

        ptCloud2 = pointCloud(vp_oct); ptCloud2.Intensity = Bp_oct(:);
        figure(14); clf; NF = 11;
        pcshow(ptCloud2); set(gcf,'color','w'); set(gca, 'FontSize',NF','Color','w','XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]); colorbar; %colormap('gray');
    end
    % [strain,zps] = phase_retr_acc_mat_unwrap(zbar, f0, BU_CC, DBU_CC, BL_CC, DBL_CC, dphi);
    % [strain,zps] =  UL_Tay_rec(zbar, f0, BU_CC, DBU_CC, BL_CC, DBL_CC, Wk);
    % clf; plot(zps,strain); hold on; plot(zp,-strain_rig)

    if not(exist("borders"))
        borders = BORD(1,:);
    end
    for ib = 2:(Nalpha_layers+1)
        indb =  find((zp2 >= borders(ib-1)) & (zp2 <= borders(ib)));
        str_line(indb) = -ALPHA(2,ib-1);
    end
else

    [strain, zp] = conventional_reconstruct(zbar, dzphi.', fit_range, 0, true);

    figure(1); clf; s1 = -alpha_max - 1e-3; s2 = alpha_min +  1e-3; ix0 =  ceil(Nscanx*0.7);
    subplot(2,2,1); imagesc(x_scan_vec*1e6,zp.*1e6,-strain.'); colorbar; caxis([s1,s2]);
    subplot(2,2,2); imagesc(x_scan_vec*1e6,zp.*1e6,-strain_rig.');  colorbar; caxis([s1,s2]);
    subplot(2,2,3); plot(zp.*1e6,-strain_rig(ix0,:)); hold on; plot(zp.*1e6,-strain(ix0,:),'r--'); hold on; plot(zp2.*1e6,str_line,'k.'); colorbar; ylim([s1,s2]);
    subplot(2,2,4); plot(zp.*1e6,-mean(strain_rig)); hold on; plot(zp.*1e6,-mean(strain),'r--'); hold on; plot(zp2.*1e6,str_line,'k.');  colorbar; ylim([s1,s2]);

    err_n(dphi,dphi_rig)
end

