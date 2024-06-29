function [Ik,varargout] = eq_sig_rig(Ik_raw,f_raw,varargin)
f_raw = f_raw(:);
Nk = length(f_raw); Nx = size(Ik_raw,2); f_min = f_raw(1);  f_max = f_raw(Nk); Df = f_max-f_min;
Nyq_z = Nk/Df/4;
Nkk = Nk; %Nk*10+1;

if not(size(Ik_raw,1) == Nk)
    error('the first dimension of the field must equal the number of frequencies');
end
Nz = 5e3+1;  zmin = -Nyq_z; zmax = Nyq_z;
if (nargin >= 3)
    Nz = varargin{1};
    if nargin >= 4
        zmin = varargin{2}; zmax = varargin{3};
    else
    end
end

f_vec = linspace(f_min,f_max,Nk); f_vec = f_vec(:); 
[A_tmp,zb2] = gen_fft2(-1,f_vec,ones(Nk,1),2*[zmin,zmax],Nz); zbar = zb2./2;
% zbar = linspace(zmin,zmax,Nz); zb2 = 2.*zbar;

A = ff_tr(Ik_raw,f_raw,zb2,-1); 
Ik = gen_fft2(+1,zb2,A,f_vec); Ik = real(Ik);

perc = 0.02;
ifv1 = find((f_vec < f_min + Df*perc)); ifv2 = find((f_vec > f_max - Df*perc)); Nfv1 = length(ifv1);
ifr1 = find((f_raw < f_vec(ifv1(Nfv1))+Df*perc)); ifr2 = find((f_raw > f_vec(ifv2(1))-Df*perc));  Nfr2 = length(ifr2); Nfv2 = length(ifv2);

Ik_tmp1 = interp1(f_raw(ifr1),Ik_raw(ifr1,:,:,:,:,:),f_vec(ifv1));
Ik_tmp2 = interp1(f_raw(ifr2),Ik_raw(ifr2,:,:,:,:,:),f_vec(ifv2));

if f_vec(Nkk) > f_raw(Nk)
    Ik_Nk = spline(f_raw(ifr2),Ik_raw(ifr2,:,:,:,:,:),f_vec(Nkk));
    Ik_tmp2(Nfv2) = Ik_Nk;
end
Ik(ifv1,:,:,:,:,:) = Ik_tmp1;
Ik(ifv2,:,:,:,:,:) = Ik_tmp2;

%     clf; plot(f_vec,Ik); hold on; plot(f_raw,Ik_raw,'--');
varargout{1} = f_vec;
varargout{2} = A;
varargout{3} = zbar;
end