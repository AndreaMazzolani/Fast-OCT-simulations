
function [EJ0_T,EJ2_T] = Ex_Debye_freq_rad_der_x_z(vertices,t0,f_vec,refind,W_0,P,K0,Th,t_w,Z_hat,varargin)

Th = Th(:); t_w = t_w(:);
LS = length(Th);

N10 = nargin == 10;
N11 = nargin == 11;
N12 = nargin == 12;

if N10 
    Nderz = 0; Nderx = 0;
elseif N11
    Nderz = varargin{1}; Nderx = 0;
elseif N12 
    Nderz = varargin{1}; Nderx = varargin{2};
else
    error(['The number of input arguments is  wrong!']);
end

c = 2.997924580105029e+08; Nk = length(f_vec); f2 = 36e-3;
w1 = pi*W_0/c;

X = vertices(:,1).';
Y = vertices(:,2).';
Z = vertices(:,3).';
Nvert = length(X);

costh = cos(Th); sinth = sin(Th);
S2 = sinth.^2;
P1 = 1 + sqrt(1-S2); P2 = (sqrt(1-S2)-1);
Jac = abs(sinth.*costh);


Rho_0 = sqrt(X.^2+Y.^2);
Phi_0 = atan2(Y,X);

sth = 2.*pi.*refind./c.*sinth;
STR = sth.*Rho_0;
Mem1 = whos('STR').bytes;
Mem = Mem1.*(1+Nderx).*(1+Nderz);

desktop = '/home/andrea/MATLAB_26_04_2022/Matlab_drive_scripts/'; % Desktop personal folder
remote = '/home/amazzolani/Matlab_drive_scripts/'; % personal folder in Mnemosine/zeus

if isfolder(desktop); M0 = 3e9;
elseif isfolder(remote); M0 = 10e9; 
else; M0 = 5e8;
end

Nm = ceil(Mem/M0);
indV = round(linspace(1,Nvert,Nm+1));
EJ0_T = zeros(Nk,Nvert,Nderz+1,Nderx+1); EJ2_T = EJ0_T;


for im = 1:Nm
    if Nm > 1;  disp(['Iter ', num2str(im),'/',num2str(Nm)]); end;
    IV = indV(im):indV(im+1); LV = length(IV);
    for ik = 1:Nk
        freq = f_vec(ik); ifc = 1i.*2.*pi.*refind.*freq./c;
        if isempty(t0) & isempty(P) 
            Mod = 1;
            if not(isempty(K0));
                cond_k0 = 1; lambda0 = 2*pi/K0; freq0 = c/lambda0;
            else 
                cond_k0 = 0;
            end
        else
            lambda0 = 2*pi/K0; freq0 = c/lambda0;
            Mod = 1i.*P.*exp(1i.*2*pi*freq*t0).*exp(-pi*P^2*(freq-freq0).^2);
        end

        XS = STR(:,IV).*freq;
        J0 = besselj(0,XS); %J0_T(ik,IV,1) = J0;
        J1 = besselj(1,XS);
        iX0 = XS == 0;
        J0_T = zeros(LS,LV,Nderx); J2_T = zeros(LS,LV,Nderx);
        J2 = 2.*1./XS.*J1 - J0; J2(iX0) = 0;

        STH = reshape(sth.*freq,[LS,1]);
        J0_T(:,:,1) = J0; J2_T(:,:,1) = J2;


        for iderx = 0:Nderx+2
            eval(['J', num2str(iderx), ' = besselj(iderx,XS);']);
        end

 
        for iderx = 1:Nderx
%             eval(['J',num2str(iderx+2),' = 2.*',num2str(iderx+1),'./XS.*J', num2str(iderx+1),' - J', num2str(iderx),';']);
%             eval(['J',num2str(iderx+2),'(iX0) = 0;']);
            id2 = ceil(iderx/2)-1; % The last term for J0 in the first loop for (The last for J2 is id2 -1;
            DJn = zeros(LS,LV); DJn2 = DJn;

            for idk = 0:(id2-1)
                DJn   = DJn   + (-1).^(idk+(iderx-2*idk)).*eval(['J', num2str(iderx-2*idk)]).*nchoosek( iderx , idk );
                DJn2  = DJn2  + (-1).^(idk+(iderx-2*idk-2)).*eval(['J', num2str(iderx-2*idk-2)]).*nchoosek( iderx , idk );
            end
            idk = id2; DJn   = DJn   + (-1).^(idk+(iderx-2*idk)).*eval(['J', num2str(iderx-2*idk)]).*nchoosek( iderx , idk );
            idk = id2; DJn2  = DJn2  + (-1).^(idk).*eval(['J', num2str(2*idk+2-iderx)]).*nchoosek( iderx , idk );
            for idk = (id2+1):iderx
                DJn  = DJn  + (-1).^(idk).*eval(['J', num2str(2*idk-iderx  )]).*nchoosek( iderx , idk );
                DJn2 = DJn2 + (-1).^(idk).*eval(['J', num2str(2*idk+2-iderx)]).*nchoosek( iderx , idk );
            end
%             DJn = DJn.*(STH.^iderx).*2.^(-iderx); DJn2 = DJn2.*(STH.^iderx).*2.^(-iderx);
            DJn = exp(log(DJn)+iderx.*log(STH/2)); DJn2 = exp(log(DJn2)+iderx.*log(STH/2)); 
            J0_T(:,:,iderx+1) = DJn; J2_T(:,:,iderx+1) = DJn2;
        end


      for iderz = 0:Nderz
          if cond_k0
              Fw = exp(-(w1.*freq0).^2.*S2).*(((costh-1).*ifc).^iderz).*exp(ifc.*(Z(IV)-Z_hat).*(costh))./(1-S2).^(1/4);
          else
            Fw = exp(-(w1.*freq).^2.*S2).*(((costh-1).*ifc).^iderz).*exp(ifc.*(Z(IV)-Z_hat).*(costh))./(1-S2).^(1/4);
          end
            TermJac = 2.*pi.*Mod.*(-1i.*f2.*freq./2./c).*Jac;
            Exs0 = TermJac.*(P1.*J0); Exs2 = TermJac.*(-P2.*J2);

            Ex0_tmp = 2*(t_w.')*Exs0./2; Ex2_tmp = 2*(t_w.')*Exs2./2;
            EJ0_T(ik,IV,iderz+1,1) = Ex0_tmp; EJ2_T(ik,IV,iderz+1,1) = Ex2_tmp;

            for iderx = 0:Nderx
                DJn = squeeze(J0_T(:,:,iderx+1)); DJn2 = squeeze(J2_T(:,:,iderx+1));
                Exs0 = TermJac.*(P1.*DJn).*Fw; Exs2 = TermJac.*(-P2.*DJn2).*Fw;

                Ex0_tmp = 2*(t_w.')*Exs0./2; Ex2_tmp = 2*(t_w.')*Exs2./2;
                EJ0_T(ik,IV,iderz+1,iderx+1) = Ex0_tmp; EJ2_T(ik,IV,iderz+1,iderx+1) = Ex2_tmp;
            end
      end

        %     dx = 7e-6; DX = dx.*ones(size(XS)); J0p = besselj(0,XS+DX.*STH); J2p = besselj(2,XS+DX.*STH);
        %     TDX = reshape((dx.^(0:Nder)./factorial(0:Nder)),[1,1,Nder+1]);
        %     J00 = sum(J0_T.*TDX,3);  J22 = sum(J2_T.*TDX,3);
        %     [err_n(J00,J0p)
        %     err_n(J0,J0p)
        %     err_n(J22,J2p)
        %     err_n(J2,J2p)
        % ]
    end

end


disp(['']);

end

