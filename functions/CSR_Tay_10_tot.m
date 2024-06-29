function [alpha_est, varargout] = CSR_Tay_10_tot(zmin,zmax,f_vec ,IU, IL, varargin) % AU, DAU, AL, DAL, Wk,  varargin)
% This function esimates the strain from IU and IL, by calculating the
% lateral displacement in "x" and "y" directions in a sparse grid, first and
% later it uses this information to calculate the strain by using all A-scans.

sIk = size(IU); NIk = length(sIk);  if NIk <= 2; Ny = 1; else; Ny = sIk(3); end;
Nx = sIk(2);

fit_range = 100e-6; res = 5e-6; 4e-6; Nx_lat = 0; Ny_lat = 0; Nx_s = 0; Ny_s = 0; Nx_s_sp = round(Nx/15); Ny_s_sp = round(Ny/15); res_sp = 20e-6;

Max_nargin = 14;
if (nargin >= 6) & (nargin <= Max_nargin); fit_range = varargin{1}; end
if (nargin >= 7) & (nargin <= Max_nargin); res = varargin{2}; end
if (nargin >= 8) & (nargin <= Max_nargin); Nx_lat = varargin{3}; end
if (nargin >= 9) & (nargin <= Max_nargin); Ny_lat = varargin{4}; end
if (nargin >= 10) & (nargin <= Max_nargin); Nx_s = varargin{5}; end
if (nargin >= 11) & (nargin <= Max_nargin); Ny_s = varargin{6}; end
if (nargin >= 12) & (nargin <= Max_nargin); Nx_s_sp = varargin{7}; end
if (nargin >= 13) & (nargin <= Max_nargin); Ny_s_sp = varargin{8}; end
if (nargin >= 14) & (nargin <= Max_nargin); res_sp = varargin{9}; end
if not( (nargin >= 5) & (nargin<= Max_nargin)); error('Wrong number of input data!'); end

Nx2_lat = Nx_lat*2; Ny2_lat = Ny_lat*2;

if abs(mod(Nx_lat,1)) + abs(mod(Ny_lat,1)) > 0; error('The lateral shifting must be an integer number!'); end

% I_norm = (mean(abs(IU)) + mean(abs(IL)))./2; IU = IU./I_norm;  IL = IL./I_norm;

f0 = mean(f_vec); f_max = max(f_vec);
cut_off = 1e-2; Wk = 2*pi*(f_max-f0)/sqrt(log(1/cut_off));
Nk = length(f_vec); Df = f_vec(Nk)-f_vec(1);

d_Disp_prec =  (sqrt(3)/2)./2./f0/11;  Disp_min_prec = -0.3e-6; Disp_max_prec = 0.3e-6; NDisp = ceil((Disp_max_prec-Disp_min_prec)./d_Disp_prec/2)*2+1;
d_alpha_prec =  0.7e-4;  alpha_min_prec = -0.3e-3; alpha_max_prec = 0.3e-3; Nalpha = ceil((alpha_max_prec-alpha_min_prec)./d_alpha_prec/2)*2+1;
vec_Disp_prec = linspace(Disp_min_prec,Disp_max_prec,NDisp); vec_alpha_prec =  linspace(alpha_min_prec,alpha_max_prec,Nalpha)+sqrt(11)*1e-9;

MagA = [1,9,33]; MagD = [1,9,33]; Niter = 3;

if not((Nx > Nx2_lat) & (Ny >= Ny2_lat)); error(['The lateral shifting in x an/or y direction it larger than the related dimension! ']); end

Ns = length(sIk); if ((Ns == 2) & (sIk(Ns) == 1)); sIk(Ns) = []; Ns = Ns-1; end % I remove the unuseful dimension
zmax_s = zmax + res_sp + 1e-6;

dz_tmp = 7e-7*sqrt(2); Nz = ceil((zmax_s - zmin)/dz_tmp/2)*2+1;
[A_tmp,zbar2] = gen_fft2(-1,f_vec(:),IU(:,1,1,1,1,1),2*[zmin ,zmax_s ],Nz);  zbar = zbar2(:)./2; dz = zbar(2) - zbar(1);

f_vec = f_vec(:);  P = 1./2./Wk.^2;
Sk =  1/sqrt(pi)./Wk.*exp(-(2*pi*(f_vec-f0)./Wk).^2);
AU = gen_fft2(-1,f_vec,IU.*Sk,zbar2); AL = gen_fft2(-1,f_vec,IL.*Sk,zbar2);
DAU = gen_fft2(-1,f_vec,(1i.*4.*pi.*f_vec).*IU.*Sk,zbar2); DAL = gen_fft2(-1,f_vec,(1i.*4.*pi.*f_vec).*IL.*Sk,zbar2);
D2AU = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^2).*IU.*Sk,zbar2); D2AL = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^2).*IL.*Sk,zbar2);
D3AU = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^3).*IU.*Sk,zbar2); D3AL = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^3).*IL.*Sk,zbar2);
D4AU = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^4).*IU.*Sk,zbar2); D4AL = gen_fft2(-1,f_vec,((1i.*4.*pi.*f_vec).^4).*IL.*Sk,zbar2);


costA = -1i.*4.*pi.*f0; Ph_zbar = exp(costA.*zbar);
BU = Ph_zbar.*AU;
DBU = Ph_zbar.*(costA.*AU + DAU);
D2BU = Ph_zbar.*((costA.^2).*AU + 2.*costA.*DAU + D2AU);
D3BU = Ph_zbar.*((costA.^3).*AU + 3.*(costA.^2).*DAU + 3.*costA.*D2AU + D3AU);
D4BU = Ph_zbar.*((costA.^4).*AU + 4.*(costA.^3).*DAU + 6.*(costA.^2).*D2AU + 4.*costA.*D3AU  + D4AU);

BL = Ph_zbar.*AL;
DBL = Ph_zbar.*(costA.*AL + DAL);
D2BL = Ph_zbar.*((costA.^2).*AL + 2.*costA.*DAL + D2AL);
D3BL = Ph_zbar.*((costA.^3).*AL + 3.*(costA.^2).*DAL + 3.*costA.*D2AL + D3AL);
D4BL = Ph_zbar.*((costA.^4).*AL + 4.*(costA.^3).*DAL + 6.*(costA.^2).*D2AL + 4.*costA.*D3AL  + D4AL);
Nb = ceil(fit_range./dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      Estimation of the lateral displacements for each point of the sparse grid

skip_sp = ceil(res_sp/dz); fit_num = min(ceil(fit_range/dz),Nz-skip_sp+1);

[Int_tmp,ind_p_sp] = slit_mean(rand(Nz,1),fit_num,1,skip_sp);
% fprintf('Nz = %d, fit_num = %d, skip_sp = %d\n',Nz,fit_num,skip_sp)
Np_sp = length(ind_p_sp); zp_sp = zbar(ind_p_sp);
vec_x_sp = (Nx_lat+1):(Nx_s_sp+1):(Nx-Nx_lat); Nx_red_sp = length(vec_x_sp);
vec_y_sp = (Ny_lat+1):(Ny_s_sp+1):(Ny-Ny_lat); Ny_red_sp = length(vec_y_sp);


% global time_thin_grid; tic;
Nx_lat0 = 0; Nx2_lat0 = 2*Nx_lat0; Nx_red0 = Nx- 2*Nx_lat0;
Ny_lat0 = 0; Ny2_lat0 = 2*Ny_lat0; Ny_red0 = Ny- 2*Ny_lat0;

Nz = Nz-skip_sp;
skip = ceil(res/dz); fit_num = min(ceil(fit_range/dz),Nz-skip+1);

[Int_tmp,ind_p] = slit_mean(rand(Nz,1),fit_num,1,skip);

Np = length(ind_p); zp = zbar(ind_p);
vec_x_T = (Nx_lat+1):(Nx-Nx_lat); Nx_red_T = length(vec_x_T);
vec_y_T = (Ny_lat+1):(Ny-Ny_lat); Ny_red_T = length(vec_y_T);

vec_x = (Nx_lat+1):(Nx_s+1):(Nx-Nx_lat); Nx_red = length(vec_x);
vec_y = (Ny_lat+1):(Ny_s+1):(Ny-Ny_lat); Ny_red = length(vec_y);

if Nx_lat+Ny_lat>0;

    % [Ib,Ialpha,IDisp,Ix_red,Iy_red,Ix_lat,Iy_lat] = ndgrid(1:Nb,1:Nalpha, 1:NDisp, 1:Nx_red_sp, 1:Ny_red_sp,1:(Nx2_lat+1),1:(Ny2_lat+1));
    Err_opt_lat_sp = zeros(Np_sp,Nx_red_sp,Ny_red_sp); Lat_Disp_x_sp = Err_opt_lat_sp; Lat_Disp_y_sp = Err_opt_lat_sp;
    alpha_est_sp = Err_opt_lat_sp; Disp_est_sp = Err_opt_lat_sp;

    Err_mat_lat_sp = zeros(1,Nx_red_sp,Ny_red_sp,Nx2_lat+1,Ny2_lat+1);
    vec_lat_x = ones(size(Err_mat_lat_sp)).*reshape(1:Nx2_lat+1,[1,1,1,Nx2_lat+1,1]);
    vec_lat_y = ones(size(Err_mat_lat_sp)).*reshape(1:Ny2_lat+1,[1,1,1,1,Ny2_lat+1]);

    ones_aD = ones(1,Nalpha,NDisp,Nx_red_sp,Ny_red_sp,Nx2_lat+1,Ny2_lat+1);
    Mat_alpha_prec = reshape(vec_alpha_prec,[1,Nalpha,1,1,1,1,1]).*ones_aD;
    Mat_Disp_prec = reshape(vec_Disp_prec,[1,1,NDisp,1,1,1,1]).*ones_aD;

    for ip = 1:Np_sp
        fprintf( 'Sparse grid: Lateral estimation: Depth index (%d/%d)\n',ip,Np_sp );

        izb_min = min(find(zbar >= zp_sp(ip) - fit_range./2)); izb_max = izb_min + Nb - 1; IZb = izb_min:izb_max; zb = zbar(IZb); zb2 = zb.*2;
        I0 = zeros(Nb,1,1,Nx_red_sp,Ny_red_sp,Nx2_lat+1,Ny2_lat+1);
        BUpre = I0; BLpre = I0; DBUpre = I0; DBLpre = I0; D2BUpre = I0; D2BLpre = I0; D3BUpre = I0; D3BLpre = I0; D4BUpre = I0; D4BLpre = I0;

        for ix0 = 1:(Nx2_lat+1) % Indices of dataset translation in x direction
            ix = ix0 - 1 - Nx_lat; shift_x = ix + vec_x_sp;
            for iy0 = 1:Ny2_lat+1
                iy = iy0 - 1 - Ny_lat; shift_y = iy + vec_y_sp;
                BUpre(:,1,1,:,:,ix0,iy0)    =   BU(IZb,vec_x_sp,vec_y_sp);    BLpre(:,1,1,:,:,ix0,iy0)   =   BL(IZb,shift_x,shift_y);
                DBUpre(:,1,1,:,:,ix0,iy0)   =  DBU(IZb,vec_x_sp,vec_y_sp);   DBLpre(:,1,1,:,:,ix0,iy0)   =  DBL(IZb,shift_x,shift_y);
                D2BUpre(:,1,1,:,:,ix0,iy0)  = D2BU(IZb,vec_x_sp,vec_y_sp);  D2BLpre(:,1,1,:,:,ix0,iy0)   = D2BL(IZb,shift_x,shift_y);
                D3BUpre(:,1,1,:,:,ix0,iy0)  = D3BU(IZb,vec_x_sp,vec_y_sp);  D3BLpre(:,1,1,:,:,ix0,iy0)   = D3BL(IZb,shift_x,shift_y);
                D4BUpre(:,1,1,:,:,ix0,iy0)  = D4BU(IZb,vec_x_sp,vec_y_sp);  D4BLpre(:,1,1,:,:,ix0,iy0)   = D4BL(IZb,shift_x,shift_y);
            end
        end

        for iter = 1:Niter
            Mat_alpha = Mat_alpha_prec.*(MagA(Niter-iter+1)); Mat_Disp = Mat_Disp_prec.*(MagD(Niter-iter+1));

            if iter == 1;
                alphaU = Mat_alpha./(2+Mat_alpha); DispU = Mat_Disp./(2+Mat_alpha);
            else
                Mat_alpha_prec_opt = alpha_opt + Mat_alpha; alphaU = Mat_alpha_prec_opt./(2+ Mat_alpha_prec_opt);
                Mat_Disp_prec_opt = Disp_opt + Mat_Disp; DispU = Mat_Disp_prec_opt./(2+ Mat_alpha_prec_opt);
            end

            alphaL = -alphaU; DispL = -DispU; costU = 1./(1+alphaU); costL = 1./(1+alphaL);
            d0U = (alphaU.*zb+DispU).^2; d2U = alphaU.*(2+alphaU); PSIU = -4.*pi.*f0.*(DispU+alphaU.*zb); Ph_BU = exp(1i.*PSIU -d0U.*Wk.^2);
            d0L = (alphaL.*zb+DispL).^2; d2L = alphaL.*(2+alphaL); PSIL = -4.*pi.*f0.*(DispL+alphaL.*zb); Ph_BL = exp(1i.*PSIL -d0L.*Wk.^2);

            c1U = 1i.*4.*pi.*f0.*alphaU + 2.*Wk.^2.*(1+alphaU).*(alphaU.*zb+DispU);
            c2U = (-Wk.^2.*d2U + c1U.^2./2);
            c3U = (-Wk.^2.*d2U.*c1U + c1U.^3./6);
            c4U = ((-Wk.^2.*d2U).^2./2 + c1U.^2.*(-Wk.^2.*d2U)./2 + c1U.^4./24);

            c1L = 1i.*4.*pi.*f0.*alphaL + 2.*Wk.^2.*(1+alphaL).*(alphaL.*zb+DispL);
            c2L = (-Wk.^2.*d2L + c1L.^2./2);
            c3L = (-Wk.^2.*d2L.*c1L + c1L.^3./6);
            c4L = ((-Wk.^2.*d2L).^2./2 + c1L.^2.*(-Wk.^2.*d2L)./2 + c1L.^4./24);

            % Calculations on Bitpaper fwhm pag. 12 -->    https://bitpaper.io/go/fwhm_W0/SJlds0iTH
            BUp1 = -P.*DBUpre;
            BUp2 = P.*(BUpre + P.*D2BUpre);
            BUp3 = -(P.^2).*(3.*DBUpre + P.*D3BUpre);
            BUp4 = (P.^2).*(3.*BUpre + 6.*P.*D2BUpre + (P.^2).*D4BUpre);

            BLp1 = -P.*DBLpre;
            BLp2 = P.*(BLpre + P.*D2BLpre);
            BLp3 = -(P.^2).*(3.*DBLpre + P.*D3BLpre);
            BLp4 = (P.^2).*(3.*BLpre + 6.*P.*D2BLpre + (P.^2).*D4BLpre);

            BUp = Ph_BU.*(BUpre  + c1U.*BUp1   + c2U.*BUp2 + c3U.*BUp3 + c4U.*BUp4);
            BLp = Ph_BL.*(BLpre  + c1L.*BLp1   + c2L.*BLp2 + c3L.*BLp3 + c4L.*BLp4);
            A1 =  costU.*(BUp.*(costA.*(zb-DispU) -1)) + costL.*(BLp.*(costA.*(zb-DispL) -1));
            A2 = costA.*(BUp + BLp); %bup = (BLp-BUp);
            BULp = BLp-BUp;

            Int_11 = sum(real(conj(A1).*A1)); Int_12 = sum(real(conj(A1).*A2));  Int_21 = Int_12;  Int_22 = sum(real(conj(A2).*A2)); det_A = (Int_11.*Int_22 - Int_12.*Int_21);
            Int_b1 = sum(real(conj(A1).*BULp)); Int_b2 = sum(real(conj(A2).*BULp));
            Int_normb = sum(conj(BULp).*(BULp));  Int_normbbu = sum(conj(BLp).*(BLp));

            eps_alpha = ( Int_22.*Int_b1 - Int_12.*Int_b2)./det_A; eps_Disp = (-Int_21.*Int_b1 + Int_11.*Int_b2)./det_A;

            alphaU_opt =  alphaU + eps_alpha; DispU_opt = DispU + eps_Disp;
            alpha_opt_vec =  2.*alphaU_opt./(1-alphaU_opt); Disp_opt_vec =  2.*DispU_opt./(1-alphaU_opt);

            Err_mat =  sqrt(real(eps_alpha.^2.*Int_11 + eps_Disp.^2.*Int_22 + eps_alpha.*eps_Disp.*(Int_12 + Int_21) - 2.*(eps_alpha.*Int_b1 + eps_Disp.*Int_b2) + Int_normb)./Int_normbbu);
            [min_alpha,Ind_alpha_opt] = min(Err_mat,[],[2,3],'linear');

            alpha_opt = alpha_opt_vec(Ind_alpha_opt); Disp_opt = Disp_opt_vec(Ind_alpha_opt);
        end

        Err_opt = Err_mat(Ind_alpha_opt);
        [Err_opt_lat_ip,in_lat_ip] = min(Err_opt,[],[6,7],'linear');
        Err_opt_lat_sp(ip,:,:) = reshape(Err_opt_lat_ip,[1,Nx_red_sp,Ny_red_sp]);
        Lat_Disp_x_sp(ip,:,:) = reshape(vec_lat_x(in_lat_ip),[1,Nx_red_sp,Ny_red_sp]);
        Lat_Disp_y_sp(ip,:,:) = reshape(vec_lat_y(in_lat_ip),[1,Nx_red_sp,Ny_red_sp]);

        alpha_est_sp(ip,:,:) =  reshape(alpha_opt(in_lat_ip),[1,Nx_red_sp,Ny_red_sp]);
        Disp_est_sp(ip,:,:) = reshape(Disp_opt(in_lat_ip),[1,Nx_red_sp,Ny_red_sp]);
    end

    Lat_Disp_x_sp = Lat_Disp_x_sp - Nx_lat-1;
    Lat_Disp_y_sp = Lat_Disp_y_sp - Ny_lat-1;


    % Extension of the arrays  "Lat_Disp_x_sp" "Lat_Disp_y_sp"
    vx_i = vec_x_sp(1); vx_f = vec_x_sp(Nx_red_sp); vec_x_sp_ext = [vx_i-(Nx_s_sp+1),vec_x_sp,vx_f+(Nx_s_sp+1)];
    vy_i = vec_y_sp(1); vy_f = vec_y_sp(Ny_red_sp); vec_y_sp_ext = [vy_i-(Ny_s_sp+1),vec_y_sp,vy_f+(Ny_s_sp+1)];
    vz_i = zp_sp(1); vz_f = zp_sp(Np_sp); zp_sp_ext = [vz_i-dz*skip_sp;zp_sp;vz_f+dz*skip_sp];

    Lat_Disp_x_sp_ext = zeros(Np_sp+2,Nx_red_sp+2,Ny_red_sp+2); Lat_Disp_x_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_x_sp;
    Lat_Disp_y_sp_ext = zeros(Np_sp+2,Nx_red_sp+2,Ny_red_sp+2); Lat_Disp_y_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_y_sp;

    Lat_Disp_x_sp_ext(1,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_x_sp(1,:,:); Lat_Disp_x_sp_ext(Np_sp+2,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_x_sp(Np_sp,:,:);
    Lat_Disp_x_sp_ext(2:Np_sp+1,1,2:Ny_red_sp+1) = Lat_Disp_x_sp(:,1,:); Lat_Disp_x_sp_ext(2:Np_sp+1,Nx_red_sp+2,2:Ny_red_sp+1) = Lat_Disp_x_sp(:,Nx_red_sp,:);
    Lat_Disp_x_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,1) = Lat_Disp_x_sp(:,:,1); Lat_Disp_x_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,Ny_red_sp+2) = Lat_Disp_x_sp(:,:,Ny_red_sp);
    Lat_Disp_y_sp_ext(1,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_y_sp(1,:,:); Lat_Disp_y_sp_ext(Np_sp+2,2:Nx_red_sp+1,2:Ny_red_sp+1) = Lat_Disp_y_sp(Np_sp,:,:);
    Lat_Disp_y_sp_ext(2:Np_sp+1,1,2:Ny_red_sp+1) = Lat_Disp_y_sp(:,1,:); Lat_Disp_y_sp_ext(2:Np_sp+1,Nx_red_sp+2,2:Ny_red_sp+1) = Lat_Disp_y_sp(:,Nx_red_sp,:);
    Lat_Disp_y_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,1) = Lat_Disp_y_sp(:,:,1); Lat_Disp_y_sp_ext(2:Np_sp+1,2:Nx_red_sp+1,Ny_red_sp+2) = Lat_Disp_y_sp(:,:,Ny_red_sp);

    vnz = [1,Np_sp+2]; vnz_sp = [1,Np_sp];
    vnx = [1,Nx_red_sp+2]; vnx_sp = [1,Nx_red_sp];
    vny = [1,Ny_red_sp+2]; vny_sp = [1,Ny_red_sp];

    for it1 = 1:2
        for it2 = 1:2
            Lat_Disp_x_sp_ext(vnz(it1),vnx(it2),2:Ny_red_sp+1)  = Lat_Disp_x_sp(vnz_sp(it1),vnx_sp(it2),:);
            Lat_Disp_x_sp_ext(2:Np_sp+1,vnx(it1),vny(it2))  = Lat_Disp_x_sp(:,vnx_sp(it1),vny_sp(it2),:);
            Lat_Disp_x_sp_ext(vnz(it1),2:Nx_red_sp+1,vny(it2))  = Lat_Disp_x_sp(vnz_sp(it1),:,vny_sp(it2));
            Lat_Disp_y_sp_ext(vnz(it1),vnx(it2),2:Ny_red_sp+1)  = Lat_Disp_y_sp(vnz_sp(it1),vnx_sp(it2),:);
            Lat_Disp_y_sp_ext(2:Np_sp+1,vnx(it1),vny(it2))  = Lat_Disp_y_sp(:,vnx_sp(it1),vny_sp(it2),:);
            Lat_Disp_y_sp_ext(vnz(it1),2:Nx_red_sp+1,vny(it2))  = Lat_Disp_y_sp(vnz_sp(it1),:,vny_sp(it2));

            for it3 = 1:2
                Lat_Disp_x_sp_ext(vnz(it1),vnx(it2),vny(it3))  = Lat_Disp_x_sp(vnz_sp(it1),vnx_sp(it2),vny_sp(it3));
                Lat_Disp_y_sp_ext(vnz(it1),vnx(it2),vny(it3))  = Lat_Disp_y_sp(vnz_sp(it1),vnx_sp(it2),vny_sp(it3));
            end
        end
    end


    % [IZ_sp,IX_sp,IY_sp] = meshgrid(zp_sp_ext.',vec_x_sp_ext,vec_y_sp_ext);
    % [IZ_sp_r,IX_r,IY_sp_r] = meshgrid(zp',vec_x,vec_y);

    vx = round(interp3(vec_x_sp_ext',zp_sp_ext.',vec_y_sp_ext,Lat_Disp_x_sp_ext,vec_x_T',zp',vec_y_T,'linear'));
    vy = round(interp3(vec_x_sp_ext',zp_sp_ext.',vec_y_sp_ext,Lat_Disp_y_sp_ext,vec_x_T',zp',vec_y_T,'linear'));

else
    vx = zeros(Np,Nx_red_T,Ny_red_T); vy = zeros(Np,Nx_red_T,Ny_red_T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Estimation of the strain, by using the results from the calculation made on the sparse grid.

% [Ib,Ialpha,IDisp,Ix_red,Iy_red,Ix_lat,Iy_lat] = ndgrid(1:Nb,1:Nalpha, 1:NDisp, 1:Nx_red, 1:Ny_red,1:(Nx2_lat+1),1:(Ny2_lat+1));
Err_opt_lat = zeros(Np,Nx_red,Ny_red); Lat_Disp_x = Err_opt_lat; Lat_Disp_y = Err_opt_lat;
alpha_est = Err_opt_lat; Disp_est = Err_opt_lat;

Err_mat_lat = zeros(1,Nx_red,Ny_red,Nx2_lat0+1,Ny2_lat0+1);
vec_lat_x = ones(size(Err_mat_lat)).*reshape(1:Nx2_lat0+1,[1,1,1,Nx2_lat0+1,1]);
vec_lat_y = ones(size(Err_mat_lat)).*reshape(1:Ny2_lat0+1,[1,1,1,1,Ny2_lat0+1]);

ones_aD = ones(1,Nalpha,NDisp,Nx_red,Ny_red,Nx2_lat0+1,Ny2_lat0+1);
Mat_alpha_prec = reshape(vec_alpha_prec,[1,Nalpha,1,1,1,1,1]).*ones_aD;
Mat_Disp_prec = reshape(vec_Disp_prec,[1,1,NDisp,1,1,1,1]).*ones_aD;


BL_old = BL;  DBL_old = DBL;  D2BL_old = D2BL;  D3BL_old = D3BL;  D4BL_old = D4BL;
for ip = 1:Np
    izb_min = min(find(zbar >= zp(ip) - fit_range./2)); izb_max = izb_min + Nb - 1; IZb = izb_min:izb_max;
    for ix = 1:Nx_red_T
        ix0 = vec_x_T(ix);
        for iy = 1:Ny_red_T
            iy0 = vec_y_T(iy);
            ivx = vx(ip,ix,iy); ivy = vy(ip,ix,iy);
            BL(IZb,ix0,iy0) = BL_old(IZb,ix0+ivx,iy0+ivy);
            DBL(IZb,ix0,iy0) = DBL_old(IZb,ix0+ivx,iy0+ivy);
            D2BL(IZb,ix0,iy0) = D2BL_old(IZb,ix0+ ivx,iy0+ivy);
            D3BL(IZb,ix0,iy0) = D3BL_old(IZb,ix0+ivx,iy0+ivy);
            D4BL(IZb,ix0,iy0) = D4BL_old(IZb,ix0+ivx,iy0+ivy);
        end
    end
end

if not((Nx_lat0 <= Nx_lat) & (Ny_lat0 <= Ny_lat)); error('The number of the improved lateral displacement must be smaller than the "sparse grid" lateral displacement!'); end


for ip = 1:Np
    fprintf( 'Thin grid: Depth index (%d/%d)\n',ip,Np );

    izb_min = min(find(zbar >= zp(ip) - fit_range./2)); izb_max = izb_min + Nb - 1; IZb = izb_min:izb_max; zb = zbar(IZb); zb2 = zb.*2;
    I0 = zeros(Nb,1,1,Nx_red,Ny_red,Nx2_lat0+1,Ny2_lat0+1);
    BUpre = I0; BLpre = I0; DBUpre = I0; DBLpre = I0; D2BUpre = I0; D2BLpre = I0; D3BUpre = I0; D3BLpre = I0; D4BUpre = I0; D4BLpre = I0;

    for ix0 = 1:(Nx2_lat0+1) % Indices of dataset translation in x direction
        ix = ix0 - 1 - Nx_lat0; shift_x = ix + vec_x;
        for iy0 = 1:Ny2_lat0+1
            iy = iy0 - 1 - Ny_lat0; shift_y = iy + vec_y;
            BUpre(:,1,1,:,:,ix0,iy0)    =   BU(IZb,vec_x,vec_y);    BLpre(:,1,1,:,:,ix0,iy0)   =   BL(IZb,shift_x,shift_y);
            DBUpre(:,1,1,:,:,ix0,iy0)   =  DBU(IZb,vec_x,vec_y);   DBLpre(:,1,1,:,:,ix0,iy0)   =  DBL(IZb,shift_x,shift_y);
            D2BUpre(:,1,1,:,:,ix0,iy0)  = D2BU(IZb,vec_x,vec_y);  D2BLpre(:,1,1,:,:,ix0,iy0)   = D2BL(IZb,shift_x,shift_y);
            D3BUpre(:,1,1,:,:,ix0,iy0)  = D3BU(IZb,vec_x,vec_y);  D3BLpre(:,1,1,:,:,ix0,iy0)   = D3BL(IZb,shift_x,shift_y);
            D4BUpre(:,1,1,:,:,ix0,iy0)  = D4BU(IZb,vec_x,vec_y);  D4BLpre(:,1,1,:,:,ix0,iy0)   = D4BL(IZb,shift_x,shift_y);
        end
    end

    for iter = 1:Niter
        Mat_alpha = Mat_alpha_prec.*(MagA(Niter-iter+1)); Mat_Disp = Mat_Disp_prec.*(MagD(Niter-iter+1));

        if iter == 1;
            alphaU = Mat_alpha./(2+Mat_alpha); DispU = Mat_Disp./(2+Mat_alpha);
        else
            Mat_alpha_prec_opt = alpha_opt + Mat_alpha; alphaU = Mat_alpha_prec_opt./(2+ Mat_alpha_prec_opt);
            Mat_Disp_prec_opt = Disp_opt + Mat_Disp; DispU = Mat_Disp_prec_opt./(2+ Mat_alpha_prec_opt);
        end

        alphaL = -alphaU; DispL = -DispU; costU = 1./(1+alphaU); costL = 1./(1+alphaL);
        d0U = (alphaU.*zb+DispU).^2; d2U = alphaU.*(2+alphaU); PSIU = -4.*pi.*f0.*(DispU+alphaU.*zb); Ph_BU = exp(1i.*PSIU -d0U.*Wk.^2);
        d0L = (alphaL.*zb+DispL).^2; d2L = alphaL.*(2+alphaL); PSIL = -4.*pi.*f0.*(DispL+alphaL.*zb); Ph_BL = exp(1i.*PSIL -d0L.*Wk.^2);

        c1U = 1i.*4.*pi.*f0.*alphaU + 2.*Wk.^2.*(1+alphaU).*(alphaU.*zb+DispU);
        c2U = (-Wk.^2.*d2U + c1U.^2./2);
        c3U = (-Wk.^2.*d2U.*c1U + c1U.^3./6);
        c4U = ((-Wk.^2.*d2U).^2./2 + c1U.^2.*(-Wk.^2.*d2U)./2 + c1U.^4./24);

        c1L = 1i.*4.*pi.*f0.*alphaL + 2.*Wk.^2.*(1+alphaL).*(alphaL.*zb+DispL);
        c2L = (-Wk.^2.*d2L + c1L.^2./2);
        c3L = (-Wk.^2.*d2L.*c1L + c1L.^3./6);
        c4L = ((-Wk.^2.*d2L).^2./2 + c1L.^2.*(-Wk.^2.*d2L)./2 + c1L.^4./24);

        % Calculations on Bitpaper fwhm pag. 12 -->    https://bitpaper.io/go/fwhm_W0/SJlds0iTH
        BUp1 = -P.*DBUpre;
        BUp2 = P.*(BUpre + P.*D2BUpre);
        BUp3 = -(P.^2).*(3.*DBUpre + P.*D3BUpre);
        BUp4 = (P.^2).*(3.*BUpre + 6.*P.*D2BUpre + (P.^2).*D4BUpre);

        BLp1 = -P.*DBLpre;
        BLp2 = P.*(BLpre + P.*D2BLpre);
        BLp3 = -(P.^2).*(3.*DBLpre + P.*D3BLpre);
        BLp4 = (P.^2).*(3.*BLpre + 6.*P.*D2BLpre + (P.^2).*D4BLpre);

        %         fprintf('size Ph_BU = %d/n size BUpre = %d/n size c1U = %d/n size c2U = %d/n', size(Ph_BU))
%         disp(size(Ph_BU)); disp(size(BUpre)); disp(size(c1U));
%         disp(size(BUp1));disp(size(c2U));disp(size(BUp2))
        BUp = Ph_BU.*(BUpre  + c1U.*BUp1   + c2U.*BUp2 + c3U.*BUp3 + c4U.*BUp4);
        BLp = Ph_BL.*(BLpre  + c1L.*BLp1   + c2L.*BLp2 + c3L.*BLp3 + c4L.*BLp4);

        %         BUp = Ph_BU.*(BUpre  + c1U.*BUp1   ); BLp = Ph_BL.*(BLpre  + c1L.*BLp1  );
        %         BUp = Ph_BU.*(BUpre    ); BLp = Ph_BL.*(BLpre  );

        A1 =  costU.*(BUp.*(costA.*(zb-DispU) -1)) + costL.*(BLp.*(costA.*(zb-DispL) -1));
        A2 = costA.*(BUp + BLp); %bup = (BLp-BUp);
        BULp = BLp-BUp;

        Int_11 = sum(real(conj(A1).*A1)); Int_12 = sum(real(conj(A1).*A2));  Int_21 = Int_12;  Int_22 = sum(real(conj(A2).*A2)); det_A = (Int_11.*Int_22 - Int_12.*Int_21);
        Int_b1 = sum(real(conj(A1).*BULp)); Int_b2 = sum(real(conj(A2).*BULp));
        Int_normb = sum(conj(BULp).*(BULp));  Int_normbbu = sum(conj(BLp).*(BLp));

        eps_alpha = ( Int_22.*Int_b1 - Int_12.*Int_b2)./det_A; eps_Disp = (-Int_21.*Int_b1 + Int_11.*Int_b2)./det_A;

        alphaU_opt =  alphaU + eps_alpha; DispU_opt = DispU + eps_Disp;
        alpha_opt_vec =  2.*alphaU_opt./(1-alphaU_opt); Disp_opt_vec =  2.*DispU_opt./(1-alphaU_opt);

        Err_mat =  sqrt(real(eps_alpha.^2.*Int_11 + eps_Disp.^2.*Int_22 + eps_alpha.*eps_Disp.*(Int_12 + Int_21) - 2.*(eps_alpha.*Int_b1 + eps_Disp.*Int_b2) + Int_normb)./Int_normbbu);
        [min_alpha,Ind_alpha_opt] = min(Err_mat,[],[2,3],'linear');

        alpha_opt = alpha_opt_vec(Ind_alpha_opt); Disp_opt = Disp_opt_vec(Ind_alpha_opt);
    end

    Err_opt = Err_mat(Ind_alpha_opt);
    [Err_opt_lat_ip,in_lat_ip] = min(Err_opt,[],[6,7],'linear');
    Err_opt_lat(ip,:,:) = reshape(Err_opt_lat_ip,[1,Nx_red,Ny_red]);
    Lat_Disp_x(ip,:,:) = reshape(vec_lat_x(in_lat_ip),[1,Nx_red,Ny_red]);
    Lat_Disp_y(ip,:,:) = reshape(vec_lat_y(in_lat_ip),[1,Nx_red,Ny_red]);

    alpha_est(ip,:,:) =  reshape(alpha_opt(in_lat_ip),[1,Nx_red,Ny_red]);
    Disp_est(ip,:,:) = reshape(Disp_opt(in_lat_ip),[1,Nx_red,Ny_red]);
end


IXT = 1:(Nx_s+1):(Nx_red_T);
IYT = 1:(Ny_s+1):(Ny_red_T);

Lat_Disp_x = Lat_Disp_x - Nx_lat0-1 + vx(:,IXT,IYT);
Lat_Disp_y = Lat_Disp_y - Ny_lat0-1 + vy(:,IXT,IYT);
time_thin_grid = toc;
varargout{1} = zp;
varargout{2} = Disp_est;
varargout{3} = Err_opt_lat;
varargout{4} = Lat_Disp_x;
varargout{5} = Lat_Disp_y;
varargout{6} = vec_x;
varargout{7} = vec_y;

if 0

    iz1 = 1;  vec_Disp_est = vec_Disp(1) + 13./2./f0 - zp(iz1).*(vec_alpha-vec_alpha(1));

    Err1 = squeeze(Err_mat2(:,:,iz1,1,1));
    figure(1); clf;
    subplot(1,3,1); plot(vec_alpha,min(Err1,[],2));
    subplot(1,3,2); plot(vec_Disp,min(Err1,[],1));
    subplot(1,3,3);  imagesc_set(Err1.',vec_Disp,vec_alpha); hold on; plot(vec_alpha,vec_Disp_est,'k*');


    % surf(Err1,vec_alpha,vec_Disp*1e6)

    [ALPHA1,DISP1] = meshgrid(vec_alpha,vec_Disp); ALPHA = ALPHA1.'; DISP = DISP1.'; clf; surf(ALPHA,DISP,Err1)

    ind_alpha_opt = zeros(Np,1); ind_Disp_opt = zeros(Np,1);
    for ip =1:Np
        Err1 = squeeze(Err_mat2(:,:,ip,1,1));
        [min_E,ind_E] = min_array(Err1); ind_alpha_opt(ip) = ind_E(1); ind_Disp_opt(ip) = ind_E(2);
        % alpha_opt(ip) = vec_alpha(ind_alpha_opt(ip)); Disp_opt(ip) = vec_Disp(ind_Disp_opt(ip));
    end
    alpha_opt = vec_alpha(ind_alpha_opt); Disp_opt = vec_Disp(ind_Disp_opt);

    ip0 = ceil(Np*0.7); vec_Disp_all = vec_Disp  - Disp_opt(ip0);
    clf; plot(vec_Disp_all.*2*f0,squeeze(Err1([ind_alpha_opt(ip0),1],:)));
    % ia = ind_alpha_opt; iD = ind_Disp_opt;
end







end











