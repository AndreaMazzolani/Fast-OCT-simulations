function [strain_UWA,varargout] = strain_UWA(lambda,zbar,dphi,x_scan_vec,varargin)

remote = '/home/amazzolani/Matlab_drive_scripts/My_papers_code/General_compression_OCE_simulations/';
desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/My_papers_code/General_compression_OCE_simulations/';
laptop = 'C:/Users/andre/Documents/MATLAB/Matlab_drive_scripts/'; % personal folder in my Laptop (October 2022)
work = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % personal folder in my Laptop (October 2022)

if isfolder(remote)
    out_dir = remote;
elseif isfolder(desktop)
    out_dir = desktop;
elseif isfolder(laptop)
    out_dir = laptop;
elseif isfolder(work)
    out_dir = work;
end

cond_out = (nargout > 1);

if (nargin <= 3) | (nargin >= 10); error(['Wrong number of input data']); end

diffz = diff(zbar); dz = diffz(1);
y_scan_vec = 0; fit_range = 100e-6; NZU = ceil(25e-6/dz);
UNW = [4,4,NZU]; AP = [3,3,3];

if nargin >= 5; y_scan_vec = varargin{1}; end
if nargin >= 6; fit_range = varargin{2}; end
if nargin >= 7; UNW =  varargin{3}; end; if not(length(UNW)==3); error('Unwrapping params must be 3'); end
if nargin >= 8; AP =  varargin{4}; end;  if not(length(UNW)==3); error('Averaging params must be 3'); end

LD_x = UNW(1); LD_y = UNW(2); LD_z = UNW(3);
Ax = AP(1); Ay = AP(2); Az = AP(3);

% Input AP -> from [z,x,y] -> [x,y,z], or [x,y,z] -> [y,z,x]  
dphi = strain_optimized2(dphi,ones(size(dphi)),Ay,Az,Ax);

s_dphi = size(dphi);
Lx = length(x_scan_vec);
Ly = length(y_scan_vec);
Lz = length(zbar);

% fprintf('Nx = %d, Ny = %d, Nz = %d\n',Lx,Ly,Lz); disp(s_dphi);

if not((Lx == s_dphi(1))&(Ly == s_dphi(2))&(Lz == s_dphi(3)))
    fprintf('s_dphi = [%d,%d,%d]\n[Lx,Ly,Lz] = [%d,%d,%d]\n',s_dphi,[Lx,Ly,Lz]); 
    error(['Wrong input data dimensions']); 
    
end

if length(x_scan_vec) > 1
    diffx = diff(x_scan_vec); dx = diffx(1);
    if not(abs(max(diffx)-min(diffx)) < 1e-10);
        error(['x_scan_vec must be equally spaced']);
    end
else dx = 1e-6;
end

if Ly == 1
    dy = 1e-6;
else
    diffy = diff(y_scan_vec); dy = diffy(1);

    if not(abs(max(diffy)-min(diffy)) < 1e-10);
        error(['y_scan_vec must be equally spaced']);
    end

end

if fit_range >  zbar(end)- zbar(1); fit_range = zbar(end)- zbar(1); end

if not(abs(max(diffz)-min(diffz)) < 1e-10);
    error(['zbar must be equally spaced']);
end

Pd.voxel_size_metres.x = dx;
Pd.voxel_size_metres.y = dy;
Pd.voxel_size_metres.z = dz;
Pd.cplx_phase_diff_xyz = dphi;

UZ = min(LD_z,Lz); UY = min(LD_y,Ly); UX = min(LD_x,Lx);
unwrap_config.PHASE_UNWRAP_Z_RANGE  = UZ;
unwrap_config.PHASE_UNWRAP_X_RADIUS = UX;
unwrap_config.PHASE_UNWRAP_Y_RADIUS = UY;
unwrap_config.PHASE_UNWRAP_ALWAYS   = true;

% out_dir = 'UWA_Testing_dataset';
% out_dir = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/My_papers_code/General_compression_OCE_simulations/';
out_name_oce  = fullfile(out_dir, ['WPU_unwrapped_phase_difference_sajat', '.mat']);

out_pd_uw = volume_unwrap_phase_lateral( out_name_oce, Pd, unwrap_config);

S1 = sqrt(abs(dphi)); % weight of OCT SNR
phase1_flip = out_pd_uw.pd_unwrapped_xyz;
 % physical distance in air
pix_size_z = dz; % axial pixel size in air
strain_fit_pix = round( fit_range / pix_size_z );
D_i  = double( pix_size_z * (1:strain_fit_pix)' );
D_i2 = double( D_i .^ 2 );


OCT = (S1(:,:,:));
var_xz = 1./OCT(:,:,:);
Disp1 = (lambda/4/pi).*phase1_flip(:,:,:);
disp_xyz1 = Disp1;

%% Calculate the size of the output arrays and strain
out_size_y = size(disp_xyz1, 2);
out_size_x = size(disp_xyz1, 1);
out_size_z = size(disp_xyz1, 3) - strain_fit_pix + 1;
iz1 = ceil(strain_fit_pix/2);
iz2 = Lz - iz1; 
z_UWA = linspace(zbar(iz1),zbar(iz2),out_size_z);


WLSLocalStrain1 = zeros(out_size_x, out_size_y, out_size_z);
% disp_xz1 = zeros(out_size_x,512);

parfor i = 1:out_size_y
    for LatInc = 1:out_size_x
        %% For each fit segment
        % Have some extra accumulation variables to make parfor happy
        disp_xz1 = squeeze(disp_xyz1(:,i,:));

        WLS_Ascan1       = zeros(out_size_z, 1);
        WLSVarEst_Ascan1 = zeros(out_size_z, 1);

        disp_Ascan1 = double( disp_xz1(LatInc,:) );
        disp_var_Ascan1 = double( var_xz(LatInc,i,:) ); % every B-scan
        for PixInc = 1:out_size_z
            % column vector of displacements that make up the fit segment
            N_disp1 = disp_Ascan1( PixInc:PixInc+(strain_fit_pix-1) );
            % variance of the displacements that make up the fit segment
            N_DispVar1 = squeeze(disp_var_Ascan1( PixInc:PixInc+(strain_fit_pix-1)))';

            % sum of the displacement variances
            k0_wls1 = sum(1./N_DispVar1);
            k1_wls1 = sum(D_i'./N_DispVar1);
            k2_wls1 = sum(D_i2'./N_DispVar1);

            %% From Kennedy 2012b:
            % WLSnumerator: k0_wls * sum(weight .* D_i .* N_disp) - k1_wls *
            %               sum(weight .* N_disp)
            % weight = 1 ./ N_DispVar;
            % WLScorrect = k0_wls * sum(weight .* D_i .* N_disp) - k1_wls * ...
            %     sum(weight .* N_disp);

            %% Simplification from Sze Howe's Honours Thesis (Koh2011)
            % units of 1 / length ^ 2
            WLSnumerator1 = sum(((k0_wls1.*D_i' - k1_wls1) ./ N_DispVar1) .* N_disp1 );

            % Converted to length units of (Kennedy 2012b) "Strain estimation in
            % phase-sensitive optical coherence elastography"

            %% From Kennedy 2012b:
            % Units of 1 / length ^ 2
            % (\sum w_j) (\sum w_j (z_j - z_{i-1})^2) - (\sum w_j (z_j - z_{i-1}))^2
            WLSdenominator1 = k0_wls1 * k2_wls1 - k1_wls1^2;

            % Unitless
            WLS_Ascan1(PixInc) = WLSnumerator1 / WLSdenominator1;

            %% From Press2007NRTASC, section 15.2
            % Estimate the variance of the WLS estimate
            WLSVarEst_Ascan1(PixInc) = k0_wls1 / WLSdenominator1;
        end
        WLSLocalStrain1(LatInc,i,:) = WLS_Ascan1;
        %         WLSVarEst(LatInc,:)      = WLSVarEst_Ascan;
    end
%     fprintf('output B-scan %d\n',i);

end

strain_UWA = WLSLocalStrain1;
if cond_out
    varargout{1} = z_UWA;
    varargout{2} = phase1_flip;
end

end

