function [a,b,c,varargout] = gauss_regrW(x,yT,varargin)
    % Gaussian interpolation:  Solution of min{||y - a*exp(-(1/2)*((x-b)/c)^2) ||} , output = a,b,c. The regression is performed along the first dimension
    
    numvar = length(varargin); 
    if numvar > 1; error('Too many input data'); end
    
    sx = size(x); Nx = sx(1); Dx = max(x)-min(x); sx_ext = sx; sx_ext(1) = Nx+2;
    x_ext = zeros(sx_ext); x_ext(2:Nx+1,:,:,:,:,:,:,:) = x; 
    x_ext(1,:,:,:,:,:,:,:) = min(x) - 5.*Dx; x_ext(Nx+2,:,:,:,:,:,:,:) = max(x) + 5.*Dx;
    My = max(abs(yT)); y = yT./My; sy = size(yT); sy_ext = sy; sy_ext(1) = Nx+2; 
    y_ext = zeros(sy_ext); y_ext(2:Nx+1,:,:,:,:,:,:,:) = y; 
    miny = min(y_ext(:)); disp = 1e-6; ind_neg = y_ext <= disp; 
    if miny <= 0; % ind_max = y == 1; % xmax = x(ind_max); %xmax = xmax(:); 
%         fprintf('\nthe "y" variable has non-positive values, they will be removed for the gaussian regression with "%e"\n\n',disp); 
        y_ext(ind_neg) = y_ext(ind_neg).*0 + 10;
%         x_quad(ind_neg) = x_quad(ind_neg).*0 + xmax(1);
    end
    y_ext(1,:,:,:,:,:,:,:) = 1e-6; y_ext(Nx+2,:,:,:,:,:,:,:) = 1e-6; w0 = 0;
%     y_ext = y; x_ext = x;
    logy = log(y_ext); W = abs(y_ext).^1.5; if numvar == 1; W(2:Nx+1,:,:,:,:,:,:,:) = varargin{1}; end; W(ind_neg) = 0.*W(ind_neg); W(1,:,:,:,:,:,:,:) = w0; W(Nx+2,:,:,:,:,:,:,:) = w0; W = W./sum(W); 
    [a0,b0,c0,y0] = quad_interp(x_ext,logy,W);
    
    cond_NAN = c0>=0; c0(cond_NAN) = NaN; 
%     if max(c0(:)) > 0; error('There are complex parameters, check the input!!'); end
    
    c = 1./sqrt(-2.*c0); %if not(isequal(c,real(c))); error('There are complex parameters, check the input!!'); end
    b = -b0./2./c0;
    a = My.*exp(a0 - (b0.^2)./4./c0);
    num_nan = max([sum(isnan(a(:))),sum(isnan(b(:))),sum(isnan(c(:)))]);
%     if num_nan>0; fprintf('There are "%d" NaN of "%d" total values, the regression in those points failed!\n', num_nan,numel(a)); end
 
    if nargout >= 4; %A = reshape(a,[1,size(x(sx(2:end)))]); B = reshape(b,[1,size(a)]); C = reshape(c,[1,size(a)]); 
        y_est = reshape(a.*exp(-((x-b).^2./2./(c.^2))),size(yT)); varargout{2} = y_est;
        y_est(:,cond_NAN) = 1; yT(:,cond_NAN) = 0; varargout{1} = err_n(y_est,yT); 
    end
    if nargout >= 6; error('The maximum variables in output must be 5'); end

end


% TEST 

% Nx = 1001; N1 = 12; N2 = 3; N3 = 4;
% A0 = rand(1,N1,N2,N3); B0 = rand(1,N1,N2,N3); C0_tmp = rand(1,N1,N2,N3); C0 = C0_tmp + 0.1;
% X = sort((2*rand(Nx,N1,N2,N3)-1)*1); Y0 =  A0.*exp(-(1/2).*((X-B0)./C0).^2);
 
% [A,B,C,E,Y] = gauss_regrW(X,Y0);
% [err_n(A,A0),err_n(B,B0),err_n(C,C0),err_n(Y,Y0)]



