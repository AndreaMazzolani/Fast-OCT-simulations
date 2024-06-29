%TO-DO
function [a,b,c,varargout] = gauss_regrW2(x,yT,varargin)
    % Gaussian interpolation:  Solution of min{||y - a*exp(-(1/2)*((x-b)/c)^2) ||} , output = a,b,c. The regression is performed along the first dimension
    
    numvar = length(varargin); 
    if numvar > 1; error('Too many input data'); end
    
    sx = size(x); Nx = sx(1); Dx = max(x)-min(x); sx_ext = sx; sx_ext(1) = Nx+2;
    x_ext = zeros(sx_ext); x_ext(2:Nx+1,:,:,:,:,:,:) = x; 
    x_ext(1,:,:,:,:,:,:) = min(x) - 5.*Dx; x_ext(Nx+2,:,:,:,:,:,:) = max(x) + 5.*Dx;
    My = max(abs(yT)); y = yT./My; sy = size(yT); sy_ext = sy; sy_ext(1) = Nx+2; 
    y_ext = zeros(sy_ext); y_ext(2:Nx+1,:,:,:,:,:,:) = y; 
    miny = min(y(:)); disp = 1e-25; ind_neg = y <= disp; 
    if miny <= 0; % ind_max = y == 1; % xmax = x(ind_max); %xmax = xmax(:); 
%         fprintf('\nthe "y" variable has non-positive values, they will be removed for the gaussian regression with "%e"\n\n',disp); 
        y(ind_neg) = y(ind_neg).*0 + 1;
%         x_quad(ind_neg) = x_quad(ind_neg).*0 + xmax(1);
    end
    y_ext(1,:,:,:,:,:,:) = 1e-6; y_ext(Nx+2,:,:,:,:,:,:) = 1e-6; 
%     y_ext = y; x_ext = x;
    logy = log(y_ext); W = y_ext.^1.5; if numvar == 1; W = varargin{1};  end; W = W./sum(W); W(ind_neg) = 0.*W(ind_neg); 
    [a0,b0,c0,y0] = quad_interp(x_ext,logy,W);
    
    cond_NAN = c0>=0; c0(cond_NAN) = NaN; 
%     if max(c0(:)) > 0; error('There are complex parameters, check the input!!'); end
    
    c = 1./sqrt(-2.*c0); %if not(isequal(c,real(c))); error('There are complex parameters, check the input!!'); end
    b = -b0./2./c0;
    a = My.*exp(a0 - (b0.^2)./4./c0);
    num_nan = max([sum(isnan(a(:))),sum(isnan(b(:))),sum(isnan(c(:)))]);
    if num_nan>0; fprintf('There are "%d" NaN of "%d" total values, the regression in those points failed!\n', num_nan,numel(a)); end
 
    if nargout >= 4; %A = reshape(a,[1,size(x(sx(2:end)))]); B = reshape(b,[1,size(a)]); C = reshape(c,[1,size(a)]); 
        y_est = reshape(a.*exp(-((x-b).^2./2./(c.^2))),size(yT)); varargout{2} = y_est;
        y_est(:,cond_NAN) = 1; yT(:,cond_NAN) = 0; varargout{1} = err_n(y_est,yT); 
    end
    if nargout >= 6; error('The maximum variables in output must be 2'); end

end


% TEST 
% 
% 
% x = linspace(-4,4,1001).';
% 
% A0 = rand(1,12,3,4); B0 = rand(1,12,3,4); C0 = rand(1,12,3,4);
% y = A0.*ones(size(x)) + B0.*x + C0.*x.^2;
% 
% C1 = C0 + 0.2;
% y2 = abs(A0).*exp(-((x-B0).^2)./2./(C1.^2));
% 
% [a,b,c] = quad_interp(x,y);
% [ag,bg,cg] = gauss_interp(x,y2);
% 
% [err_n(ag,abs(A0))
% err_n(bg,B0)
% err_n(cg,C1),
% err_n(a,A0)
% err_n(b,B0)
% err_n(c,C0)]
% 


