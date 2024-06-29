function [Fy,varargout] = gen_fft2(coef,t,F,fv,varargin)

%This function calculates an accurate and fast approximation of the "Continuous" Fourier
%transform (or the inverse Fourier transform) by employing the built-in
%function "fft", for arbitrary sampled domain and codomain intervals, the
%only assumption is that they are "equally-spaced"

%DEFINITIONS:
%         Fourier transform: -->  "F(v):= Integral(f(t)*exp(-i*2*pi*t*v)*dt) "
% Inverse Fourier transform: -->  "F^(-1)(v) := Integral(f(t)*exp(+i*2*pi*t*v)*dt) "


%INPUT:
%          coef -> ("+1","-1"); It specifies if the function will calculate the Fourier transform ("coef = +1") or its inverse(" coef = -1").
%            t  -> Unidimensional array of equally spaced sampled points of the domain variable.
%            F  -> Sampled function of t. (The first dimension of F must equal the dimension of "t")
%            fv ->  -(if fv is a number -> fv = f0= "CENTRAL FREQUENCY". In this case the interval of frequencies is automatically calculated in the function)
%                   -(if fv is an interval of points, the function will evaluate the Fourier transform in the interval "fv"
%                   -(if fv is made by two numbers (fv = [fv1,fv2]), the
%                   function will evaluate the Fourier transform in the
%                   interval "linspace(fv1,fv2,Nf)", where Nf is the number
%                   of frequencies
% "varargin{1}" -> It is the number of frequencies(Nf) where the Fourier transform is evaluated. It is taken in account by the function only if "fv" is composed by one number or 2 numbers. In that case, "varargin{1} = number  of frequencies where the Fourier transform will be
%                  evaluated. If it is missed -> "varargin{1} = length(t)".

%OUTPUT:
%             Fy ->  Fourier transform( or inverse Fourier transform) of the sampled array F.
% "varargout{1}" ->  Interval of points where the Fourier transform (or its inverse) is evaluated.
%                    If "fv" is an array ->  "varargout{1} = fv",  if "fv" is a number -> "varargout{1}" is an interval of points centered in fv.

%NOTE: The fastest and most accurate case are when "fv" is a number or two numbers, which represents the central frequency
%      where the Fourier transform(or its inverse) will be evaluated.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = 1e7; N_acc = 5e4; %Maximum number of frequencies we let calculate by "fft" Matlab function (More this number is large and more the approximation is accurate but slower)
floor_eps = 1e-6; L2_cond = 0;


if (max(abs(diff(diff(t))))>1e-8*max(abs(diff(t))))
    error('The domain has not been equally sampled');
end
t = t(:); Nt = length(t);
if isempty(varargin)
    Nf = Nt;
elseif length(varargin) == 1
    Nf = varargin{1};
elseif length(varargin) == 2
    if ((class(varargin{1}) == 'char') & strcmp(varargin{1},'pad'))
    else
        error('The 5-th input must be the string "pad" in the case of 6 input')
    end

else
    error('There are too many input arguments');
end
if ((coef ~= 1) & (coef ~= -1))
    error('"coef"  must be "+1" or "-1", which means that the function will calculate the Fourier transform or its inverse');
end

sF = size(F);
if sF(1) == length(t) % Checking that the first dimension of F equals t
    st = [0:(Nt-1)]'; t0 = mean(t); dt = abs(t(2)-t(1)); Mt = (Nt-1)/2;
    if ((length(fv) == 1)& (Nf>=Nt)) % In this case fv is a number, and the interval of frequencies is calculated inside the function
        f0 = fv;  Mf = (Nf-1)/2;
        f_vec = [-Mf:Mf]./Nf./dt+f0;  sf = [0:(Nf-1)]'; % f_vec is the array of frequencies where the the Fourier transform (or its inverse) is calculated
        Cf = exp(-1i.*coef.*2.*pi.*(f0-Mf./dt./Nf).*st.*dt);
        Cc = exp(-1i.*coef.*2*pi.*(t0-Mt.*dt).*(f0+(sf-Mf)./dt./Nf));
        if coef ==1
            Fy = fft(F.*Cf,Nf).*Cc.*dt; Nrr = Nf;
        else
            Fy = ifft(F.*Cf,Nf).*Nf.*Cc.*dt; Nrr = Nf;
        end
        varargout{1} = f_vec;
    else % In this case is necessary to find the "right number of zero-padding" in order to estimate the Fourier transform (or its inverse)

        if  (length(fv) == 1) % In conjunction with the previous "if", this statement implies even that "Nf< Nt",
            % and in this case the array of frequencies "fv" is calculated below.
            Mf = (Nf-1)/2;
            fv = [-Mf:Mf]./Nf./dt+fv;
        elseif  (length(fv) == 2) % The input is given by the first and last frequency
            f_min = fv(1); f_max = fv(2); L2_cond = 1;
            if (length(varargin) == 2) % it means that varargin(2) is the number of zero padding
                Np2 = varargin{2}; % zero-padding (chosen)
                fv = [f_min:1/dt/Np2:f_max];  Nf = length(fv); Mf = (Nf-1)/2;
            else
                Np2 = round((Nf-1)/dt/(fv(end)-fv(1))); % number of frequencies (chosen)
                f_max2 = f_min + 1/dt/Np2*(Nf-1);
                fv = [f_min:1/dt/Np2:f_max2];
                Nf = length(fv); Mf = (Nf-1)/2;
            end

        else % The input of frequencies fv is an array
            Nf = length(fv); Mf = (Nf-1)/2;
            if (max(abs(diff(diff(fv))))>1e-8*max(abs(diff(fv))))
                error('The codomain has not been equally sampled');
            end
        end
        f0 = mean(fv);  if (L2_cond); df = 1/dt/Np2; else; df = abs(fv(2)-fv(1)); end
        %The following rows estimate the minimum number of frequencies "N_opt"
        %between the frequencies "fv" such that "|fv(j)-fv(j+1)|< 10^(-2)|*fv(j)|"

        % Calculation of the "right" zero-padding just estimating p,q:
        % df*dt = p/q. In this case the denominator q is the right number of zero padding.
        x = df*dt;

        xr(1) = mod(x(1),1); j1 = 1; err_toll = 1; toll = 1e-6;

        %err_toll2 = min(mod(1/x(1),1),abs(1 - mod(1/x(1),1))); % this will check if " df ≈ 1/Nf/dt" for some N

        q1 = round(1/x(1));  err_toll2 = abs(df-1/q1/dt)/df; q = q1;

        if ( (err_toll2 < 1e-8) & (q<Nmax) ) % in that case df ≈ 1/q/dt
            p1 = 1;

            if (q < Nt) % In this case I find a suitable multiple of "q" that is greater than Nt
                qq = Nt-mod(Nt,q)+q;  %% qq = q*Nq (multiple)
                p1 = round(p1/q*qq+floor_eps); %p1 s.t. p1/q1 = p/q
                q = round(qq +floor_eps);
            end

            Nf1 = q; Mf1 = (Nf1-1)/2; %f01 = f0 + 1/dt*Mf1/Nf1- Mf*df;
            tmp = 0; NftM = 0; Fy = zeros([Nf,sF(2:end)]);

            while (tmp == 0)
                Nftm = NftM+1;
                f01 = f0 + (Mf1)/dt/Nf1 - (Mf-NftM)/dt/q1; %fvec1 = [-Mf1:Mf1]/dt/Nf1 + f01;
                sf1 = [0:(Nf1-1)]';
                Cf1 = exp(-1i.*coef.*2.*pi.*(f01-Mf1./dt./Nf1).*st.*dt);
                Cc1 = exp(-1i.*coef.*2*pi.*(t0-Mt.*dt).*(f01+(sf1-Mf1)./dt./Nf1));

                if coef == 1
                    Fy1 = fft(F.*Cf1,Nf1).*Cc1.*dt; Nrr = Nf1; %Zero padding
                else
                    Fy1 = ifft(F.*Cf1,Nf1).*Nf1.*Cc1.*dt;  Nrr = Nf1;%Zero padding
                end

                Nft1 = floor((Nf1-1)/p1 + Nftm + floor_eps); NftM = Nft1;

                if (NftM >= Nf)
                    NftM = Nf; tmp = 1;
                end

                if (NftM == Nf)& (Nftm == 1)
                    if (p1  == 1)
                            Fy = Fy1(1:NftM-Nftm+1,:,:,:,:); 
                    else
                            Fy = Fy1(p1*(0:(NftM-Nftm))+1,:,:,:,:); 
                    end
                else
                        Fy(Nftm:NftM,:,:,:,:) = Fy1(p1*(0:(NftM-Nftm))+1,:,:,:,:); 

%                     for j1 = Nftm:NftM
%                         jp = p1*(j1-Nftm)+1;
% 
%                         if length(sF) ==1
%                             Fy(j1) = Fy1(jp);
%                         elseif length(sF) ==2
%                             Fy(j1,:) = Fy1(jp,:);
%                         elseif length(sF) ==3
%                             Fy(j1,:,:) = Fy1(jp,:,:);
%                         elseif length(sF) ==4
%                             Fy(j1,:,:,:) = Fy1(jp,:,:,:);
%                         end
%                     end
                end

            end
            varargout{1} = fv;

        else
            disp(['Frequencies in gen_fft2 are calculated completely']);
            while ((err_toll(j1) > 1e-6) & (abs(xr(j1)) > 1e-7) & (j1 < 20)); %|(j1<3))
                x(j1+1) = 1/xr(j1);
                xr(j1+1) = mod(x(j1+1),1);
                j1 = j1+1;

                if (x(j1) == 0)
                    j4 = j1-1;
                else
                    j4 = j1;
                end

                rat_x(j1) = floor(x(j1) + floor_eps);

                for j2 = 2:j4
                    j3 = j4-j2+1;
                    rat_x(j1) = floor(x(j3)+floor_eps)+ 1/rat_x(j1);
                end
                err_toll(j1) = abs(rat_x(j1)-abs(x(1)))/abs(x(1));
            end

            Lrx = length(rat_x);
            if(Lrx >= 2)
                [p,q] = numden(sym(rat_x)); PP0 = double(p); QQ0 = double(q);
            else
                [p,q] = numden(sym(rat_x)); PP0 = double(p); QQ0 = double(q);
            end

            p = [double(p)]; q = [double(q)]; Np = length(p);

            %  ???           if (coef == 1) %(coef = 1 -> Fourier transform)
            %                 Ft = fft(F(:,1)); Ft2 = fft(t(:).*F);
            %             else %(coef = -1 -> Inverse Fourier transform)
            %                 Ft = ifft(F); Ft2 = ifft(t(:).*F);
            %             end
            %
            %             % Calculation of the maximum number of frequencies "N_opt" that
            %             % will be used to estimate the Fourier transform
            %             cond = (abs(Ft2) > 1e-2.*max(abs(Ft2)));
            %             TT = abs(Ft./Ft2); TT = TT(cond); h_opt = 10^(-3)/2/pi*min(TT,[],'all','linear');
            %             N_opt = round((1/dt)/h_opt+floor_eps)+1;
            %             if (N_opt>Nmax)
            %                 N_opt = Nmax;
            %             end

            %             if(max(abs(TT),[],'all') == 0 )
            %                 error('The function  is zero everywhere');
            %             end

            N_opt = Nmax;


            if Np>2 % In that case I look for a better "2°" rational approximation

                p_opt = double(p(Np)); q_opt = double(q(Np));
                p1 = double(p(Np-1)); q1 = double(q(Np-1));
            else
                NN = ceil(N_acc/q(Np));
                J2 = 0:1e3;
                if (p(Np)/q(Np) < x(1))

                    rat_x3 = (NN.*p(Np) + J2)./(NN.*q(Np));
                    Ix3 = max(find(rat_x3 < x(1)));
                    p_opt = (NN*p(Np) + J2(Ix3)); q_opt = NN*q(Np);

                elseif (p(Np)/q(Np) == x(1))
                    p_opt = p(Np); q_opt = q(Np);

                else
                    rat_x3 = (NN.*p(Np) - J2)./(NN.*q(Np));
                    Ix3 = max(find(rat_x3 > x(1)));
                    p_opt = (NN*p(Np) - J2(Ix3)); q_opt = NN*q(Np);
                end
            end

            if (q(Np) < N_opt)
                p1 = double(p(Np-1)); q1 = double(q(Np-1));
                p2 = double(p(Np)); q2 = double(q(Np));

                if p2 == 0
                    error('The frequencies are too close each other')
                end

                if p1 == 0
                    p1 = p2;
                    q1 = q2;
                end

                if (q1 < Nt)
                    qq1 = Nt - mod(Nt,q1) + q1;
                    p1 = round(p1/q1*qq1 + floor_eps);
                    q1 = round(qq1 + floor_eps);
                end

                if (q2 < Nt)
                    qq2 = Nt - mod(Nt,q2) + q2;
                    p2 = round(p2/q2*qq2 + floor_eps);
                    q2 = round(qq2 + floor_eps);
                end

                Nf1 = q1; Mf1 = (Nf1-1)/2; %f01 = f0 + 1/dt*Mf1/Nf1- Mf*df;
                Nf2 = q2; Mf2 = (Nf2-1)/2; %f02 = f0 + 1/dt*Mf2/Nf2- Mf*df;
                tmp = 0; NftM = 0;  Fy = zeros([Nf,sF(2:end)]);

                while (tmp == 0)
                    Nftm = NftM+1;
                    f01 = f0 + (Mf1)/dt/Nf1- (Mf-NftM)*df;
                    f02 = f0 + (Mf2)/dt/Nf2- (Mf-NftM)*df;
                    fvec1 = [-Mf1:Mf1]./Nf1./dt+f01;  sf1 = [0:(Nf1-1)]';
                    Cf1 = exp(-1i.*coef.*2.*pi.*(f01-Mf1./dt./Nf1).*st.*dt);
                    Cc1 = exp(-1i.*coef.*2*pi.*(t0-Mt.*dt).*(f01+(sf1-Mf1)./dt./Nf1));

                    if coef == 1
                        Fy1 = fft(F.*Cf1,Nf1).*Cc1.*dt; Nrr = Nf1;  %Zero padding
                    else
                        Fy1 = ifft(F.*Cf1,Nf1).*Nf1.*Cc1.*dt; Nrr = Nf; %Zero padding
                    end

                    fvec2 = [-Mf2:Mf2]./Nf2./dt+f02;  sf2 = [0:(Nf2-1)]';
                    Cf2 = exp(-1i.*coef.*2.*pi.*(f02-Mf2./dt./Nf2).*st.*dt);
                    Cc2 = exp(-1i.*coef.*2*pi.*(t0-Mt.*dt).*(f02+(sf2-Mf2)./dt./Nf2));

                    if coef ==1
                        Fy2 = fft(F.*Cf2,Nf2).*Cc2.*dt; Nrr = Nf2;  %Zero padding
                    else
                        Fy2 = ifft(F.*Cf2,Nf2).*Nf2.*Cc2.*dt;  Nrr = Nf2; %Zero padding
                    end

                    Nft1 = floor((Nf1-1)/p1+Nftm+floor_eps);
                    Nft2 = floor((Nf2-1)/p2+Nftm+floor_eps);

                    NftM = min(Nft1,Nft2);
                    if (NftM >= Nf)
                        NftM = Nf; tmp = 1;
                    end

                    for j1 = Nftm:NftM
                        jp1 = p1*(j1-Nftm)+1;
                        jp2 = p2*(j1-Nftm)+1;

                        if fvec2(jp2) == fvec1(jp1)
                            if length(sF) ==1
                                Fy(j1) = Fy2(jp2);
                            elseif length(sF) ==2
                                Fy(j1,:) = Fy2(jp2,:);
                            elseif length(sF) ==3
                                Fy(j1,:,:) = Fy2(jp2,:,:);
                            elseif length(sF) ==4
                                Fy(j1,:,:,:) = Fy2(jp2,:,:,:);
                            end
                            Fy(j1) = Fy2(jp2);
                        else
                            if length(sF) == 1
                                Fy(j1) = (Fy2(jp2).*(fv(j1)-fvec1(jp1))+Fy1(jp1).*(fvec2(jp2)-fv(j1)))/(fvec2(jp2)-fvec1(jp1));
                            elseif length(sF) == 2
                                Fy(j1,:) = (Fy2(jp2,:).*(fv(j1)-fvec1(jp1))+Fy1(jp1,:).*(fvec2(jp2)-fv(j1)))/(fvec2(jp2)-fvec1(jp1));
                            elseif length(sF) == 3
                                Fy(j1,:,:) = (Fy2(jp2,:,:).*(fv(j1)-fvec1(jp1))+Fy1(jp1,:,:).*(fvec2(jp2)-fv(j1)))/(fvec2(jp2)-fvec1(jp1));
                            elseif length(sF) == 4
                                Fy(j1,:,:,:) = (Fy2(jp2,:,:,:).*(fv(j1)-fvec1(jp1))+Fy1(jp1,:,:,:).*(fvec2(jp2)-fv(j1)))/(fvec2(jp2)-fvec1(jp1));
                            end
                        end
                    end
                end
                varargout{1} = fv;
            else


                p = PP0; q = QQ0;

                iq = min(find(err_toll < 3e-3));
                if isempty(iq)
                    iq = max(find(q <= N_acc));
                end


                Nf3 = double(q(iq));  Mf3 = (Nf3-1)/2;
                NftM = 0; tmp = 0; Fy = zeros([Nf,sF(2:end)]);
                while (tmp == 0)
                    Nftm = NftM+1;
                    f03 = f0 + (Mf3)/dt/Nf3- (Mf-NftM)*df;
                    fvec3 = [-Mf3:Mf3]./Nf3./dt + f03;  sf3 = [0:(Nf3-1)]';
                    Cf3 = exp(-1i.*coef.*2.*pi.*(f03-Mf3./dt./Nf3).*st.*dt);
                    Cc3 = exp(-1i.*coef.*2*pi.*(t0-Mt.*dt).*(f03+(sf3-Mf3)./dt./Nf3));

                    if coef ==1
                        Fy3 = fft(F.*Cf3,Nf3).*Cc3.*dt; Nrr = Nf3; %Zero padding
                    else
                        Fy3 = ifft(F.*Cf3,Nf3).*Nf3.*Cc3.*dt; Nrr = Nf3; %Zero padding
                    end

                    NftM = floor((Nf3-1)/(df*dt*Nf3)+Nftm+floor_eps);

                    if (NftM >= Nf)
                        NftM = Nf; tmp = 1;
                    end

                    for j1 = Nftm:NftM
                        jp1 = floor(df*dt*Nf3*(j1-Nftm)+floor_eps)+1;
                        jp2 = jp1+1;

                        if     (length(sF) == 1)
                            Fy(j1) = (Fy3(jp2).*(fv(j1)-fvec3(jp1))+Fy3(jp1).*(fvec3(jp2)-fv(j1)))/(fvec3(jp2)-fvec3(jp1));
                        elseif (length(sF) == 2)
                            Fy(j1,:) = (Fy3(jp2,:).*(fv(j1)-fvec3(jp1))+Fy3(jp1,:).*(fvec3(jp2)-fv(j1)))/(fvec3(jp2)-fvec3(jp1));
                        elseif (length(sF) == 3)
                            Fy(j1,:,:) = (Fy3(jp2,:,:).*(fv(j1)-fvec3(jp1))+Fy3(jp1,:,:).*(fvec3(jp2)-fv(j1)))/(fvec3(jp2)-fvec3(jp1));
                        elseif l(length(sF) == 4)
                            Fy(j1,:,:,:) = (Fy3(jp2,:,:,:).*(fv(j1)-fvec3(jp1))+Fy3(jp1,:,:,:).*(fvec3(jp2)-fv(j1)))/(fvec3(jp2)-fvec3(jp1));
                        end

                    end
                end
                varargout{1} = fv;
            end
        end
    end
else
    error('The array of the domain and the sample function have different size');
end
end

