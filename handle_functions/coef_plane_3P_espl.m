function [A,B,C] = coef_plane_3P_espl(Pz1,Pz2,Pz3)

% The result are the coefficients A,B,C of the plane function "z = Ax+By+C" passing for the points Pz_j

s1 = size(Pz1); s2 = size(Pz2); s3 = size(Pz3); 
if not(isequal(s1,s2) & isequal(s1,s3)); error(['The dimension of the 3 points must match;']); end
if length(s1) <= 2 & min(s1) == 1; Pz1 = Pz1(:); Pz2 = Pz2(:); Pz3 = Pz3(:); s1 = s1(:); s2 = s2(:); s3 = s3(:); end
 
if not(s1(1) == 3); error('The first dimension must be --> (x,y,z) coordinates'); end 

if isequal(Pz1,Pz2) | isequal(Pz1,Pz3); error('The three points must be different and no collinear.'); end

% P1 = Pz1(1:2,:,:,:,:,:,:); P2 = Pz2(1:2,:,:,:,:,:,:); P3 = Pz3(1:2,:,:,:,:,:,:);

V13 = Pz3-Pz1; x13 = V13(1,:,:,:,:,:,:); y13 = V13(2,:,:,:,:,:,:); z13 = V13(3,:,:,:,:,:,:);
V23 = Pz3-Pz2; x23 = V23(1,:,:,:,:,:,:); y23 = V23(2,:,:,:,:,:,:); z23 = V23(3,:,:,:,:,:,:);
V12 = Pz2-Pz1; x12 = V12(1,:,:,:,:,:,:); y12 = V12(2,:,:,:,:,:,:); z12 = V12(3,:,:,:,:,:,:);

Den = y13.*x12 - y12.*x13;

cond_parallel = Den == 0;

if sum(cond_parallel(:)) > 0; error('The 3 points on the (x,y) plane are collinear!!'); end

B = (x12.*z13-x13.*z12)./Den;
A = (z12.*y13-z13.*y12)./Den;
C = Pz1(3,:,:,:,:,:,:) - A.*Pz1(1,:,:,:,:,:,:)  - B.*Pz1(2,:,:,:,:,:,:); 



end


% TEST with a general example;
% 
% N1 = 10; N2 = 20; N3 = 4; N4 = 7;
% 
% NT = 1e3; 
% A = NT*(2*rand(1,N1,N2,N3,N4)-1);
% B = NT*(2*rand(1,N1,N2,N3,N4)-1);
% C = NT*(2*rand(1,N1,N2,N3,N4)-1);
% 
% P1  = NT*(2*rand(3,N1,N2,N3,N4)-1); % points
% coef_par = 0.7;
% 
% P2 = P1; P2(1,:,:,:,:) = P2(1,:,:,:,:)+NT.*(1+ coef_par*(2*rand(1,N1,N2,N3,N4)-1)); % I move the x coordinates of P2 respect with P1;
% Q2 = P1; Q2(2,:,:,:,:) = Q2(2,:,:,:,:)+NT.*(1+ coef_par*(2*rand(1,N1,N2,N3,N4)-1)); % I move the x coordinates of P2 respect with P1;
% Q1 = P1; Q1(2,:,:,:,:) = Q2(2,:,:,:,:)+NT.*(0+ coef_par*(2*rand(1,N1,N2,N3,N4)-1)); % I move the x coordinates of P2 respect with P1;
%          Q1(1,:,:,:,:) = P2(1,:,:,:,:)+NT.*(0+ coef_par*(2*rand(1,N1,N2,N3,N4)-1)); % I move the x coordinates of P2 respect with P1;
% P1(3,:,:,:,:) = A.*P1(1,:,:,:,:) + B.*P1(2,:,:,:,:) + C;  
% P2(3,:,:,:,:) = A.*P2(1,:,:,:,:) + B.*P2(2,:,:,:,:) + C;  
% Q1(3,:,:,:,:) = A.*Q1(1,:,:,:,:) + B.*Q1(2,:,:,:,:) + C;  
% Q2(3,:,:,:,:) = A.*Q2(1,:,:,:,:) + B.*Q2(2,:,:,:,:) + C;  
% 
% cond = (round(rand(N1,N2,N3,N4)) == 1);
% Ps1 = P1; Ps2 = P2; Qs1 = Q1; Qs2 = Q2;
% Ps1(:,cond) = Q2(:,cond); Ps2(:,cond) = Q1(:,cond); Qs1(:,cond) = P1(:,cond);  Qs2(:,cond) = P2(:,cond);
% 
% [A_est,B_est,C_est] = coef_plane_3P_espl(P1,Q1,P2); 
% m1 = max([err_n(A_est,A), err_n(B_est,B) ,err_n(C_est,C)]);
% 
% [A_est,B_est,C_est] = coef_plane_3P_espl(Ps1,Qs2,Ps2); 
% m2 = max([err_n(A_est,A), err_n(B_est,B) ,err_n(C_est,C)]);
% 
% [A_est,B_est,C_est] = coef_plane_3P_espl(Ps1,Ps2,Qs1); 
% m3 = max([err_n(A_est,A), err_n(B_est,B) ,err_n(C_est,C)]);
% 
% [m1,m2,m3]


