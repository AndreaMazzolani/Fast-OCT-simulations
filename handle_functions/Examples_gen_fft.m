% Examples of application of the function "gen_fft2". This function
% calculates an accurate approximation of the continuous Fourier transform (or its
% inverse) for general sampled domain and codomain arrays.

% There are 2 assumptions:
% 1) The domain and codomain have been sampled evenly.
% 2) The domain variable has been sampled such that the Fourier transform can be approximated by the discrete Fourier transform. 
%         Nnamely, the function approaches to zero at the edges of the sampled domain and the points in the domain encode the slope of the function.

close all; clc;  
Nt = 200; t = linspace(-5,8,Nt)'; 
dt = diff(t(1:2)); %Sampling spacing;


%Calculaion of the "NATURAL" Discrete frequencies Domain 
f_vec = (0:(Nt-1))'/Nt/dt;
f_vec = fftshift(f_vec);
iz0 = find( f_vec == 0);
f_vec( 1:(iz0-1) ) = f_vec( 1:(iz0-1) ) - f_vec(iz0-1) - f_vec(iz0+1);

% Calculation of an "arbitrary" Discrete frequencies Domain 
f_vec2 = linspace(f_vec(1),f_vec(end),455)*1.5 + mean(f_vec)*0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Discrete frequencies domain could be even calculated by :
% -->    Mt =(Nt-1)/2; df = 1/Nt/dt; f_vec = df*linspace(-Mt,Mt,Nt)'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Handle - functions of the function "iF" and its analytical Fourier transform ("F")
iF = @(tvec) exp(-pi.*((tvec-1).^2)).*(tvec-2);
F = @(freq) exp(-pi.*(freq).^2).*(-1i.*freq - 1).*exp(-1i*2*pi*freq);

%Sampled function in the interval "t"
ft = iF(t);
%% Some examples of application of the function gen_fft2

%%5 examples of Fourier transform approximation 
[Fy1,fv1]  = gen_fft2(1,t,ft,mean(f_vec));       % 1): Fourier transform of "ft" evaluated in an interval of the codomain, which is centered at the point "mean(f_vec)". In this case, the interval is set by the function according to "fft" Matlab function requirements. [dt --> 1/N/dt]
[Fy2,fv2]  = gen_fft2(1,t,ft,mean(f_vec),10*Nt); % 2): Same case of 1) but we have chosen the number of output points where the Fourier transform is evaluated (N = 10*Nt)
[Fy3,fv3]  = gen_fft2(1,t,ft,mean(f_vec),30);    % 3): Same case of 1) and 2) but we calculate the Fourier transform only in 30 points  (N = 30).
[Fy4,fv4]  = gen_fft2(1,t,ft,f_vec);             % 4): The most general case, we have chosen the interval of sampled values of the codomain where the Fourier transform is evaluated
[Fy5,fv5]  = gen_fft2(1,t,ft,f_vec2);            % 5): Same case as 4), but in a different set of points of the codomain (f_vec2)

%Rigorous Fourier transform calculated for the same cases, calculated in the intervals "fv1","fv2","fv3","fv4","fv5".
Fr1 =  F(fv1(:));
Fr2 =  F(fv2(:));
Fr3 =  F(fv3(:));
Fr4 =  F(fv4(:));
Fr5 =  F(fv5(:));

%Indices where the functions are not (practically) zero% 
I1 = (find(abs(Fy1)>1e-6)); 
I2 = (find(abs(Fy2)>1e-6)); 
I3 = (find(abs(Fy3)>1e-6)); 
I4 = (find(abs(Fy4)>1e-6)); 
I5 = (find(abs(Fy5)>1e-6));

figure(1); clf; sgtitle('Fourier transform');
subplot(1,3,1); plot(fv2(I2),real(Fr2(I2))); hold on; plot(fv1(I1),real(Fy1(I1))); hold on; plot(fv2(I2),real(Fy2(I2))); hold on; plot(fv3(I3),real(Fy3(I3)),'*'); hold on; plot(fv4(I4),real(Fy4(I4))); hold on; plot(fv5(I5),real(Fy5(I5)),'.');  title('real(Fy)'); legend('anal','1)','2)','3)','4)','5)');
subplot(1,3,2); plot(fv2(I2), abs(Fr2(I2))); hold on; plot(fv1(I1), abs(Fy1(I1))); hold on; plot(fv2(I2), abs(Fy2(I2))); hold on; plot(fv3(I3), abs(Fy3(I3)),'*'); hold on; plot(fv4(I4), abs(Fy4(I4))); hold on; plot(fv5(I5), abs(Fy5(I5)),'.');  title('real(Fy)'); legend('anal','1)','2)','3)','4)','5)');
subplot(1,3,3); plot(fv1(I1),abs(Fr1(I1)-Fy1(I1))./abs(Fr1(I1))); hold on;  plot(fv2(I2),abs(Fr2(I2)-Fy2(I2))./abs(Fr2(I2))); hold on; plot(fv3(I3),abs(Fr3(I3)-Fy3(I3))./abs(Fr3(I3))); hold on; plot(fv4(I4),abs(Fr4(I4)-Fy4(I4))./abs(Fr4(I4))); hold on; plot(fv5(I5),abs(Fr5(I5)-Fy5(I5))./abs(Fr5(I5)));  title('Relative error');  legend('1)','2)','3)','4)','5)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Inverse Fourier transform   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we calculate the inverse Fourier transform to re-calculate ft from Fy, by using ph_ifft

t2 = linspace (-2,5,543); % An arbitary  interval where to calculate the inverse Fourier transform.
[ft1,tv1]  = gen_fft2(-1,fv1,Fy1,mean(t));       % 1) Inverse Fourier transform of "Fy1" in an interval of points (in the codomain) centered at "mean(t)". In this case, the "output interval" is set automatically by the function according to "ifft" Matlab function requirements. [dt --> 1/N/dt]
[ft2,tv2]  = gen_fft2(-1,fv2,Fy2,mean(t),10*Nt); % 2): Same case of 1) but we have chosen the number of output points of the codomain (N = 10*Nt)
[ft3,tv3]  = gen_fft2(-1,fv1,Fy1,mean(t),30);    % 3): Same case of 1) and 2) but we calculate the inverse Fourier transform only in 30 points  (N = 30).
[ft4,tv4]  = gen_fft2(-1,fv4,Fy4,t);             % 4): The most general case, we have chosen the sampled values of the codomain where the function will be calculated ( in this case is "t")
[ft5,tv5]  = gen_fft2(-1,fv5,Fy5,t2);            % 5): Same case as 4), but in a different sampled values of the codomain (t2)

%Rigorous Inverse Fourier transform calculated for the same cases, calculated in the intervals "tv1","tv2","tv3","tv4","tv5".
fr1 =  iF(tv1(:));
fr2 =  iF(tv2(:));
fr3 =  iF(tv3(:));
fr4 =  iF(tv4(:));
fr5 =  iF(tv5(:));

%Indices where the functions are not (practically) zero
J1 = (find(abs(ft1)>1e-6)); 
J2 = (find(abs(ft2)>1e-6)); 
J3 = (find(abs(ft3)>1e-6)); 
J4 = (find(abs(ft4)>1e-6)); 
J5 = (find(abs(ft5)>1e-6));

figure(2); clf; sgtitle('Inverse Fourier transform');
subplot(1,3,1); plot(tv2(J2),real(fr2(J2))); hold on; plot(tv1(J1),real(ft1(J1))); hold on; plot(tv2(J2),real(ft2(J2))); hold on; plot(tv3(J3),real(ft3(J3)),'*'); hold on; plot(tv4(J4),real(ft4(J4))); hold on; plot(tv5(J5),real(ft5(J5)));  title('real(Fy)'); legend('anal','1)','2)','3)','4)','5)');
subplot(1,3,2); plot(tv2(J2), abs(fr2(J2))); hold on; plot(tv1(J1), abs(ft1(J1))); hold on; plot(tv2(J2), abs(ft2(J2))); hold on; plot(tv3(J3), abs(ft3(J3)),'*'); hold on; plot(tv4(J4), abs(ft4(J4))); hold on; plot(tv5(J5), abs(ft5(J5)));  title('real(Fy)'); legend('anal','1)','2)','3)','4)','5)');
subplot(1,3,3); plot(tv1(J1),abs(fr1(J1)-ft1(J1))./abs(fr1(J1))); hold on;  plot(tv2(J2),abs(fr2(J2)-ft2(J2))./abs(fr2(J2))); hold on; plot(tv3(J3),abs(fr3(J3)-ft3(J3))./abs(fr3(J3))); hold on; plot(tv4(J4),abs(fr4(J4)-ft4(J4))./abs(fr4(J4))); hold on; plot(tv5(J5),abs(fr5(J5)-ft5(J5))./abs(fr5(J5)));  title('Relative error');  legend('1)','2)','3)','4)','5)');




