
function [ExB,varargout] = ExB_Debye_freq_rad(vertices,t0,f_vec,refind,W_0,P,K0,Th,t_w,varargin);
Th = Th(:); t_w = t_w(:);
Nk = length(f_vec);

N10 = (nargin == 10); N11 = (nargin == 11); N09 = (nargin == 9);

if N09; Nx = 1; Ny = 1; x_scan_vec = 0; y_scan_vec = 0; end
if N10; x_scan_vec = varargin{1}; Nx = length(x_scan_vec); Ny = 1;  y_scan_vec = 0; end
if N11; x_scan_vec = varargin{1}; y_scan_vec = varargin{2}; Nx = length(x_scan_vec); Ny = length(x_scan_vec);  end
if nargin > 11; error('Too many input arguments'); end


c = 2.997924580105029e+08; f2 = 36e-3;
w1 = pi*W_0/c;

x = vertices(:,1).';
y = vertices(:,2).';
z = vertices(:,3).';
Nv = length(x);
Nv_tot = Nv*Nx*Ny;

vert_tot = zeros(Nv_tot,3);
it = 0;

for iy = 1:Ny
    for ix = 1:Nx
        it = it+1;
        vert_tot(1+(it-1)*Nv:it*Nv,1) = x - x_scan_vec(ix);
        vert_tot(1+(it-1)*Nv:it*Nv,2) = y - y_scan_vec(iy);
        vert_tot(1+(it-1)*Nv:it*Nv,3) = z;
    end
end

X = vert_tot(:,1).';
Y = vert_tot(:,2).';
Z = vert_tot(:,3).';



S2 = sin(Th).^2;
P1 = 1 + sqrt(1-S2); P2 = (sqrt(1-S2)-1);
Ex = zeros(Nk,Nv_tot); Jac = abs(sin(Th).*cos(Th));
iv = (nargout == 2);

if iv
    Ey = Ex;
end


Rho_0 = sqrt(X.^2+Y.^2);
Phi_0 = atan2(Y,X);

% STR = 2.*pi.*refind./c.*sin(Th).*Rho_0; Mem = whos('STR').bytes;
Mem = Nv_tot*length(S2)*8;

if isunix; M0 = 2e9 ; else; M0 = 5e8; end;
Nm = ceil(Mem/M0);
indV = round(linspace(1,Nv_tot,Nm+1)); 
cond_k0 = 0;

for im = 1:Nm
%     if Nm > 1;  disp(['Iter ', num2str(im),'/',num2str(Nm)]); end;
    IV = indV(im):indV(im+1);
    for ik = 1:Nk
%         disp([num2str(ik),'/',num2str(Nk)]);

        freq = f_vec(ik);
        if isempty(t0) & isempty(P)
            Mod = 1;
            if not(isempty(K0));
                lambda0 = 2*pi/K0; freq0 = c/lambda0;
            end
        else
            lambda0 = 2*pi/K0; freq0 = c/lambda0;
            Mod = 1i.*P.*exp(1i.*2*pi*freq*t0).*exp(-pi*P^2*(freq-freq0).^2);
        end

        %       XS = STR(:,IV).*freq;
        XS = 2.*pi.*refind./c.*sin(Th).*Rho_0(1,IV).*freq;
        J0 = besselj(0,XS); J2 = besselj(2,XS);

        if cond_k0
            Fw = exp(-(w1.*freq0).^2.*S2).*exp(1i.*2.*pi.*refind.*Z(IV).*freq./c.*sqrt(1-sin(Th).^2))./(1-S2).^(1/4);
        else
            Fw = exp(-(w1.*freq).^2.*S2).*exp(1i.*2.*pi.*refind.*Z(IV).*freq./c.*sqrt(1-sin(Th).^2))./(1-S2).^(1/4);
        end

        ExS = 2.*pi.*Mod.*(-1i.*f2.*freq./2./c).*Fw.*(P1.*J0-P2.*J2.*cos(2.*Phi_0(IV))).*Jac;
        Ex_tmp = 2.*(t_w.')*ExS./2; Ex(ik,IV) = Ex_tmp;

        if iv
            EyS = 2.*pi.*Mod.*(-1i.*f2.*freq./2./c).*Fw.*(-P2.*J2.*sin(2.*Phi_0(IV))).*Jac;
            Ey_tmp = 2.*(t_w.')*EyS./2; Ey(ik,IV) = Ey_tmp;
        end

        if iv
            varargout{1} = Ey;
        end
    end
end

ExB = zeros(Nk,Nv,Nx,Ny);
it = 0;
for iy = 1:Ny
    for ix = 1:Nx
        it = it+1;
        ExB(:,:,ix,iy) = Ex(:,1+(it-1)*Nv:it*Nv,1);
    end
end
ExB = squeeze(ExB);

if iv
    it = 0;
    for iy = 1:Ny
        for ix = 1:Nx
            it = it+1;
            EyB(:,:,ix,iy) = Ey(:,1+(it-1)*Nv:it*Nv,1);
        end
    end
    ExB = squeeze(ExB); varargout{1} = EyB;
end
% disp(['']);

%[mean(mean(Tp)),max(max(Tp)),var(mean(Tp)),vp,wp*max(max(S2)),mean(mean(I1)),mean(mean(I3)),mean(mean(I2)),F(i2,i1)]'
end

