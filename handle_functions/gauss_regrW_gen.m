%TO-DO
function [a,b,c,d,varargout] = gauss_regrW_gen(x,yT2,varargin)
% Gaussian interpolation:  Solution of min{||y - a*exp(-(1/2)*((x-b)/c)^2) ||} , output = a,b,c. The regression is performed along the first dimension

numvar = length(varargin);
if numvar > 1; error('Too many input data'); end
minY2 = min(yT2); maxY2 = max(yT2); init_D = 1e-12.*(maxY2-minY2) - minY2;
yT = yT2 + init_D;
sx = size(x); Nx = sx(1); Dx = max(x)-min(x); sx_ext = sx; sx_ext(1) = Nx+2;
x_ext = zeros(sx_ext); x_ext(2:Nx+1,:,:,:,:,:,:) = x; 
x_ext(1,:,:,:,:,:,:) = min(x) - 5.*Dx; x_ext(Nx+2,:,:,:,:,:,:) = max(x) + 5.*Dx;
My = max(abs(yT)); disp = 1e-4.*My; y_in = yT; sy = size(yT); sy_ext = sy; sy_ext(1) = Nx+2; Nsy = length(sy);
% y_ext = zeros(sy_ext); y_ext(2:Nx+1,:,:,:,:,:,:) = y_in;
% D_old = zeros([1,sy(2:Nsy)])+min(y_ext)-10*disp; D_new = D_old-disp;  
D_new = zeros([1,sy(2:Nsy)]); cond_GD = 1; min_err = 5e-3;

% y_old = y_in-D_old; miny = min(y_old(:));  ind_neg = y_old <= disp;
% if miny <= 0;  y_old(ind_neg) = y_old(ind_neg).*0 + disp; end; y_ext = zeros(sy_ext); y_ext(2:Nx+1,:,:,:,:,:,:) = y_old;
% y_ext(1,:,:,:,:,:,:) = disp/10; y_ext(Nx+2,:,:,:,:,:,:) = disp/10;
% logy = log(y_ext); W = y_ext.^1.5; W(ind_neg) = 0.*W(ind_neg); W(1,:,:,:,:,:,:,:) = w0;  W(Nx+2,:,:,:,:,:,:,:) = w0; if numvar == 1; W = varargin{1};  end; W = W./sum(W); 
% [a0,b0,c0,y0] = quad_interp(x_ext,logy,W); cond_NAN = c0>=0; c0(cond_NAN) = NaN;
% c = 1./sqrt(-2.*c0);  b = -b0./2./c0; a = My.*exp(a0 - (b0.^2)./4./c0);
% num_nan = max([sum(isnan(a(:))),sum(isnan(b(:))),sum(isnan(c(:)))]);
% if num_nan>0; fprintf('There are "%d" NaN of "%d" total values, the regression in those points failed!\n', num_nan,numel(a)); end
% y_est = reshape(a.*exp(-((x-b).^2./2./(c.^2)))+D_old,size(yT)); 
% err_tmp_old = 2.*sqrt(sum(abs(y_est-yT).^2))./(sqrt(sum(abs(y_est).^2))+sqrt(sum(abs(yT).^2)));  
% err_tmp_old(isnan(err_tmp_old)) = rand;
iter = 0; N_iter = 1000; Learn_err = 7e-0;
w0 = 0e-5; %D_new = D_new.*0 + 0.25;
while cond_GD
    iter = iter +1;
    y_new = y_in-D_new; miny = min(y_new(:)); ind_neg = y_new <= disp;
    if miny <= 0; error(' y should not be negative'); end  
%     if miny <= 0;  y_new(ind_neg) = y_new(ind_neg).*0 + disp(ind_neg); end; 
    y_ext = zeros(sy_ext); y_ext(2:Nx+1,:,:,:,:,:,:) = y_new;
    y_ext(1,:,:,:,:,:,:) = disp/10; y_ext(Nx+2,:,:,:,:,:,:) = disp/10;
    logy = log(y_ext); W = abs(y_ext).^1.5; W(ind_neg) = 0.*W(ind_neg); W(1,:,:,:,:,:,:,:) = w0; W(Nx+2,:,:,:,:,:,:,:) = w0; 
    if numvar == 1; W(2:Nx+1,:,:,:,:,:) = varargin{1};  end; W = W./sum(W); 
    [a0,b0,c0,logy0] = quad_interp(x_ext,logy,W); cond_NAN = c0>=0; c0(cond_NAN) = NaN;
    c = 1./sqrt(-2.*c0);  b = -b0./2./c0; a = exp(a0 - (b0.^2)./4./c0);
    num_nan = max([sum(isnan(a(:))),sum(isnan(b(:))),sum(isnan(c(:)))]);
    if num_nan>0; fprintf('There are "%d" NaN of "%d" total values, the regression in those points failed!\n', num_nan,numel(a)); end
    y_est = reshape(a.*exp(-((x-b).^2./2./(c.^2)))+D_new,size(yT));
%     err_tmp = 2.*sqrt(sum(abs(y_est-yT).^2))./(sqrt(sum(abs(y_est).^2))+sqrt(sum(abs(yT).^2)));
    err_tmp = sum((y_est-yT).^2); grad_err = 2.*sum(y_est-yT);
    err_test = 2.*sqrt(sum(abs(y_est-yT).^2))./(sqrt(sum(abs(y_est).^2))+sqrt(sum(abs(yT).^2)));
    cond_err_min = err_test > min_err;  
    
%     D_new = D_new - Learn_err.*grad_err./(1e-10+abs(grad_err));  err_tmp_old = err_tmp; sum(err_tmp(:))
    D_new = D_new - Learn_err.*grad_err; mean(err_test(:))
%       plot(x,yT,'b'); hold on; plot(x,y_est,'r'); title(num2str(D_new));
    if (iter > N_iter) | (sum(cond_err_min(:)) == 0); cond_GD = 0; end
end

if nargout >= 5; %A = reshape(a,[1,size(x(sx(2:end)))]); B = reshape(b,[1,size(a)]); C = reshape(c,[1,size(a)]);
    varargout{2} = y_est-init_D;
    y_est(:,cond_NAN) = 1; yT(:,cond_NAN) = 0; varargout{1} = err_n(y_est,yT);
end
if nargout >= 6; error('The maximum variables in output must be 2'); end

end


% TEST

% Nx = 5; N1 = 12; N2 = 3; N3 = 4;
% A0 = rand(1,N1,N2,N3); B0 = rand(1,N1,N2,N3)*3e-1; D0 = (2*rand(1,N1,N2,N3)-1); C0_tmp = rand(1,N1,N2,N3); C0 = C0_tmp + 0.2;
% X = sort((2*rand(Nx,N1,N2,N3)-1)*1); Y0 =  A0.*exp(-(1/2).*((X-B0)./C0).^2)+D0; 
% Y10 = Y1-D0;
% 
% [A,B,C,E,Y] = gauss_regrW(X,Y0-D0,abs(Y0));
% [err_n(A,A0),err_n(B,B0),err_n(C,C0),err_n(Y,Y0-D0)]
% 
% [A,B,C,E,Y] = gauss_regrW_gen(X,Y0-D0,abs(Y0));
% [err_n(A,A0),err_n(B,B0),err_n(C,C0),err_n(Y,Y0)]
% 
% y2 = y_ext(:,11,2,4); W2 = W(:,11,2,4); W2(y2==10) = []; x2= x_ext(:,11,2,4); x2(y2==10) = []; y2(y2==10) = [];
% [a0,b0,c0,y0] = quad_interp(x2,log(y2),W2);



