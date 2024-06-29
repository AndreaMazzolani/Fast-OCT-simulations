function [alpha_est, varargout] = CSR_Tay_10_tot_comb(zmin,zmax,f_vec ,IU, IL, varargin)
% (zmin,zmax,f_vec ,Sk_raw.*IU, Sk_raw.*IL,fit_range,resolution,Nx_lat,Ny_lat)
sIk = size(IU); NIk = length(sIk);  if NIk <= 2; Ny = 1; else; Ny = sIk(3); end;
Nx = sIk(2);

fit_range = 100e-6; res = 5e-6; 4e-6; Nx_lat = 0; Ny_lat = 0; Nx_s = 0; Ny_s = 0; Nx_s_sp = 5;  1; Ny_s_sp = 1;  res_sp = 20e-6;

Max_nargin = 14;
if (nargin >= 6) & (nargin <= Max_nargin); fit_range = varargin{1}; end
if (nargin >= 7) & (nargin <= Max_nargin); res = varargin{2}; end
if (nargin >= 8) & (nargin <= Max_nargin); Nx_lat = varargin{3}; end
if (nargin >= 9) & (nargin <= Max_nargin); Ny_lat = varargin{4}; end
if (nargin >= 10) & (nargin <= Max_nargin); Nx_s_sp = varargin{5}; end
if (nargin >= 11) & (nargin <= Max_nargin); Ny_s_sp = varargin{6}; end
if (nargin >= 12) & (nargin <= Max_nargin); Nx_s = varargin{7}; end
if (nargin >= 13) & (nargin <= Max_nargin); Ny_s = varargin{8}; end
if (nargin >= 14) & (nargin <= Max_nargin); res_sp = varargin{8}; end
if not( (nargin >= 5) & (nargin<= Max_nargin)); error('Wrong number of input data!'); end

 
if abs(mod(Nx_lat,1)) + abs(mod(Ny_lat,1)) > 0; error('The lateral shifting must be an integer number!'); end
if (Nx < 2*Nx_lat); error('The lateral displacement indices in x - direction are larger than the related dimension'); end
if (Ny < 2*Ny_lat); error('The lateral displacement indices in y - direction are larger than the related dimension'); end

Nx_red = Nx - 2*Nx_lat; Ny_red = Ny - 2*Ny_lat;
Nx_red_s = Nx_red/(Nx_s_sp+1); Ny_red_s  = Ny_red/(Ny_s_sp+1); 
N_TOT1 = Nx_red_s*Ny_red_s*(2*Nx_lat+1)*(2*Ny_lat+1)
N_TOT2 = Nx_red*Ny_red
N_TOT = max(N_TOT1,N_TOT2);
% Nmax = 2e4; % Clio ok (30 secs)
% Nmax = 3e4; % Clio ok (50 secs);
Nmax = 4e4; % Clio ok (64 secs);
% Nmax = 5e4; % Clio NO (> 300 secs)
% Nmax = 1e3; % laptop NO
% Nmax = 3e2; % laptop NO

Niter_d = (N_TOT/Nmax);
Niter = ceil(Niter_d);
Nx_iter = ceil(Nx/Niter_d);
cond_while = 1; it = 1;

while (it <= Niter) & cond_while
    fprintf('\n\n Iter %d/%d \n\n',it,Niter);
    IX_tmp = (it-1).*(Nx_iter) + [1:(Nx_iter+2*Nx_lat)];
    IX_tmp_red = (it-1).*(Nx_iter) + [1:Nx_iter]; %  IX_tmp(1:Nvr);
    if max(IX_tmp) > Nx;
%         if (it < Niter); error('There is an error with the subsampling!'); end
        if (it < Niter); cond_while = 0; end
        if (IX_tmp(1) > Nx-2*Nx_lat); cond_while = 0; end; %error('There is an error with the lateral indices!'); end
        IX_tmp = IX_tmp(1):Nx;
        IX_tmp_red = IX_tmp_red(1):(Nx-2*Nx_lat);
    end

    if isempty(IX_tmp_red); 
        if it == 1; error('There is an error with the arrays subsampling.'); end
        cond_while = 0;
    else
        [str_tmp , z_csr , Disp_tmp,Err_tmp,vec_lat_x_tmp,vec_lat_y_tmp] = CSR_Tay_10_tot(zmin,zmax,f_vec ,IU(:,IX_tmp,:), IL(:,IX_tmp,:),fit_range,res,Nx_lat,Ny_lat,Nx_s,Ny_s,Nx_s_sp,Ny_s_sp,res_sp);
        if it == 1;  Np = length(z_csr); alpha_est = zeros(Np,Nx_red,Ny_red); Disp_est = alpha_est; Err_opt_lat = alpha_est; Lat_Disp_x = alpha_est; Lat_Disp_y = alpha_est; end
        alpha_est(:,IX_tmp_red,:) = squeeze(str_tmp); Disp_est(:,IX_tmp_red,:) = squeeze(Disp_tmp); % alpha_est(:,IX_tmp_red,:) +1;
        Err_opt_lat(:,IX_tmp_red,:) = squeeze(Err_tmp); Lat_Disp_x(:,IX_tmp_red,:) = squeeze(vec_lat_x_tmp); Lat_Disp_y(:,IX_tmp_red,:) = squeeze(vec_lat_y_tmp);
    end
    it = it+1;

end


varargout{1} = z_csr;
varargout{2} = Disp_est;
varargout{3} = Err_opt_lat;
varargout{4} = Lat_Disp_x;
varargout{5} = Lat_Disp_y;
% varargout{6} = vec_x;
% varargout{7} = vec_y;

end

