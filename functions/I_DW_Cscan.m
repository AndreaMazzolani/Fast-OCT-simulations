% - I call A-signal ,B-signals and C-signals the frequency domain signals whose FFT gives the A-scans B-scans and C-scans
% - This function calculates the Unloaded and Loaded C-signals related to  a set of scatterers and lateral scanning of the illumination modelled by the 3D Debye-Wolf integral.
% - The scattered field is modelled assuming First-Born approximation of the scatterers

%INPUT: 
    % "v": Dimension (Nv,3). Each row is made by the 3D coordinates (x_j,y_j,z_j) of a single scatterer point. "Nv" is the number of scatterers of the simulations  
    % "rho": Dimension (Nv,1). Reflectivity coefficients of all scatterers in the simulation.
    % "f_vec_t": Dimension (1,Nk). Spectrum of frequencies of the simulation: ( f_vec_t = c/lambda_vec) . "Nk" is the number of frequencies of the simulation. 
    % "refind" = Background refractive index (It is assumed uniform)
    % "w0" = Parameter related to the waist radius of the gaussian function in the Debye-Wolf integral. It is given by  the waist radius multiplied by the ratio of the two focal lengths of the optical system (w0 = waist*f_objective/f_collimator)
    % "NA" = numerical aperture of the optical system
    % "z_f" = z-coordinate of the focus (respect to the mirror)   -->   If 'z_f > 0'  the mirror is placed at the left of the focal plane of the objective lens (in the sample space)
    % "varargin" = x_scan_vec. This is the array of the lateralignal scanning of the focussed beam, which gives the dimesion of the B-signals. If that variable is missed then "x_scan_vec = 0", and the B-signals is composed by a single A-signal.


%OUTPUT: All of the outputs are Bsignals having dimension (Nk,Nscanx), where "Nk" is the number of frequencies of the simulation, and "Nscanx" is the number of A-scans for each B-scan.
    % "I_CC" = "Cross-correlation" terms.
    % "varargout{1} = I_AC + I_CC": "Auto correlation + Cross-correlation" terms.
    % "varargout{2} = I_DC_sc": "DC terms" related to the scatterer points
    % "varargout{3} = I_DC_ref": "DC term" related to the mirror
    % "varargout{4} = ak_sc": Sample arm modal coefficient
    % "varargout{5} = ak_ref": Mirror modal coefficient
    
function [I_CC, varargout] = I_DW_Cscan(v,rho,f_vec_t,refind,w0,NA,z_f,varargin);

cond_out = (nargout > 1);
sv = size(v);
Nrho = length(rho);
Nk = length(f_vec_t);
if nargin == 7
    x_scan_vec = 0; y_scan_vec = 0;
elseif nargin == 8
    x_scan_vec = varargin{1};    
    y_scan_vec = 0;
elseif nargin == 9
    x_scan_vec = varargin{1};    
    y_scan_vec = varargin{2};    
else 
    error(['The number of input data is incorrect!']);
end

Nscanx = length(x_scan_vec); Nscany = length(y_scan_vec);
Nv = sv(1);

if not((sv(2) == 3) & (Nv == Nrho))
    error('The length of the reflectivity profiles array must equal the number of scatterers!');
end

rho = rho(:);
asap = asin(NA/refind); % Maximum angle in the Debye-Wolf integral
n_int_th = ceil((asap)./(1e-3));  % Number of points set for the Debye-Wolf numerical integral
[th,t_w] = gauss_legendre(0,asap,n_int_th);  % Numerical integral quadrature values

t0  = ''; P = ''; K0 = '';  % parameters that I need to pass to the Debye-Wolf integral function ( In this case they are empty)
rho_max = 70e-6;

ak_mirror_folder = 'ak_ref_modal_coeffs';
ak_tot_mirror_folder = [ak_mirror_folder,'/alpha_mirr_tot_',num2str(w0*1e6),'_zfr_', num2str(z_f*1e6*refind),'.mat'];
global ak_glob; 

if isempty('ak_glob')
    if exist(ak_tot_mirror_folder) == 2;
        DA = load(ak_tot_mirror_folder);
        a_ref_tot = DA.a_ref_tot; f_long_t = DA.f_long_t;
    else
        disp(['The mirror modal coefficient needs to be calculated, because this combination of  waist radius and the focal parameters is set for the first time!']);
        c = 2.997924580105029e+08; fL1 = c/300e-9; fL2 = c/2000e-9; f_long_t = linspace(fL2,fL1,3e4+1);
        L_ref = 120e-6; Nvr = round((2*L_ref)^2/(10e-6)^2*10);
        v_ref = zeros(Nvr,3); v_ref(:,1) = -L_ref + rand(Nvr,1)* 2*L_ref; v_ref(:,2) = -L_ref + rand(Nvr,1)* 2*L_ref; v_ref(:,3) = -ones(Nvr,1).*z_f;
        Ex_mirr = ExB_Debye_freq_rad(v_ref,t0,f_long_t,refind,w0,P,K0,th,t_w,0); a_ref_tot = sum(Ex_mirr.^2,2);
        save(ak_tot_mirror_folder,'a_ref_tot','f_long_t');
    end
    ak_ref0 = interp1(f_long_t,a_ref_tot,f_vec_t); % I interpolate the modal coefficient of the mirror in this interval of frequencies
    

else
    ak_ref0 = ak_glob;
end

    % I make an array 'ak_ref' having the same dimension of the sample arm modal coefficient:
ak_ref = zeros(Nk,Nscanx,Nscany);
for iscanx = 1:Nscanx
    for iscany = 1:Nscany
        ak_ref(:,iscanx,iscany) = ak_ref0(:);
    end
end

I_CC = zeros(size(ak_ref)); % Cross correlation B-scan
if cond_out
    I_DC_ref =  (ak_ref).*conj(ak_ref); % DC terms of the mirror
    I_DC_sc = zeros(size(I_DC_ref)); % DC terms of the scatterers (Inizialitazion)
    I = zeros(size(I_DC_ref));  % Entire signal: AC+CC+DC terms (Inizialitazion)
end

ak_sc = zeros(Nk,Nscanx,Nscany); % Modal coefficient of the sample arm (Inizialitazion)

x = v(:,1); y = v(:,2);
total_inds = 1:Nv;

ZF = zeros(Nv,3);
ZF(:,3) = ones(Nv,1).*z_f;
 global T_time LB NB;
 T2 = 0;
 if T_time > 0; 
     T2 = Nscanx*Nscany*(LB-1); str_UL = 'Loaded Bscan: ';
 end
str_UL = ['Cscan ', num2str(LB), '/', num2str(NB)]; 

 for iscanx = 1:Nscanx
    for iscany = 1:Nscany

        tic;
     
        % x-component of the electric field in all scatterers and all
        % frequencies, for a single lateral scanning position:
    
        v_iter = v; v_iter(:,1) = v(:,1) - x_scan_vec(iscanx); v_iter(:,2) = v(:,2) - y_scan_vec(iscany);
        rho2 = sqrt((v_iter(:,1)).^2+(v_iter(:,2)).^2); 
        ind_DW = find(rho2 <= rho_max);
    
        Ex_vec_tmp = zeros(Nk,Nrho);
        Ex_vec_tmp2 = ExB_Debye_freq_rad(v_iter(ind_DW,:)-ZF(ind_DW,:),t0,f_vec_t,refind,w0,P,K0,th,t_w); 
        Ex_vec_tmp(:,ind_DW) = Ex_vec_tmp2;
        Ex_vec = squeeze(Ex_vec_tmp).*(sqrt(rho).'); %Field renormalized with the reflectivity profiles
        Ex_vec2 = Ex_vec.^2;
        ak_sc(:,iscanx,iscany) = squeeze(sum(Ex_vec2,2)); % Modal coefficient
        I_CC(:,iscanx,iscany)  = 2.*real(ak_sc(:,iscanx,iscany).*conj(ak_ref(:,iscanx,iscany)));
    
        if cond_out
            I_DC_sc(:,iscanx,iscany) = squeeze(sum(Ex_vec2.*conj(Ex_vec2),2));
            I(:,iscanx,iscany) = (ak_sc(:,iscanx,iscany) + ak_ref(:,iscanx,iscany)).*conj(ak_sc(:,iscanx,iscany) + ak_ref(:,iscanx,iscany));
        end
    
        T_time = T_time + toc;
        iscan = iscany + Nscany*(iscanx-1);

        TOT_time = NB*Nscanx*Nscany/(iscan+T2)*(T_time);
        if T2 > 0; time_reg = 1; else; time_reg = (0.7*((Nscanx*Nscany-iscan)/Nscanx/Nscany).^4+1); end
        estim_remaining_time = (TOT_time-T_time)*time_reg;
        s_est = seconds(estim_remaining_time); 
        s_est.Format = 'hh:mm:ss';
    
        fprintf([str_UL,', A-line (',num2str(iscanx),',',num2str(iscany),')/(',num2str(Nscanx),',',num2str(Nscany),').   Estimated remaining time: ']);
        disp(string(s_est));
    end
 end
 disp(" ");
    

if cond_out
    I_DC = I_DC_sc + I_DC_ref; % DC terms
    I_AC = I - I_CC - I_DC; % Autocorreletion terms
    
    varargout{1} = I_AC + I_CC;
    varargout{2} = I_DC_sc;
    varargout{3} = I_DC_ref;
    varargout{4} = ak_sc;
    varargout{5} = ak_ref;
end

end
