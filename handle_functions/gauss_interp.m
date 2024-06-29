function [a,b,c,varargout] = gauss_interp(x,yT)
% Gaussian interpolation:  Solution of min{||y - a*exp(-(1/2)*((x-b)/c)^2) ||} , output = a,b,c. The regression is performed along the first dimension

dispy = 0;
My = max(abs(yT)); y = yT./My; if min(y(:)) <=0; disp('the "y" variable has some non-positive values!'); dispy = exp(-10)-min(y); y = y + dispy; end;

logy = log(y);

[a0,b0,c0] = quad_interp(x,logy);

Nx = length(x); NLY = length(size(y));
cond_NAN = c0>=0; 
if max(cond_NAN(:))==1; cond_NAN2 = reshape(cond_NAN,size(y,[2:NLY])); %c0(cond_NAN) = NaN;
    x_ext = zeros(Nx+2,1); x_ext(2:Nx+1) = x; x_ext(1) = -max(abs(x(:)))*5; x_ext(Nx+2) = -x_ext(1);
    logy_NAN = logy(:,cond_NAN2); logy_ext = -15.*ones(Nx+2,size(logy_NAN,2)); logy_ext(2:Nx+1,:) = logy_NAN;
    [a0_ext,b0_ext,c0_ext] = quad_interp(x_ext,logy_ext);
    a0(cond_NAN) = a0_ext; b0(cond_NAN) = b0_ext; c0(cond_NAN) = c0_ext;
end

%     if max(c0(:)) > 0; error('There are complex parameters, check the input!!'); end

c = 1./sqrt(-2.*c0); %if not(isequal(c,real(c))); error('There are complex parameters, check the input!!'); end
b = -b0./2./c0;
a = My.*(exp(a0 - (b0.^2)./4./c0)-dispy);
num_nan = max([sum(isnan(a(:))),sum(isnan(b(:))),sum(isnan(c(:)))]);
if num_nan>0; fprintf('There are "%d" NaN of "%d" total values, the regression in those points failed!\n', num_nan,numel(a)); end

if nargout >= 4; A = reshape(a,[1,size(a)]); B = reshape(b,[1,size(a)]); C = reshape(c,[1,size(a)]);
    y_est = reshape(A.*exp(-((x-B).^2./2./(C.^2))),size(y)); varargout{2} = y_est;
    y_est(:,cond_NAN) = 1; yT(:,cond_NAN) = 0; varargout{1} = err_n(y_est,yT);
    
%     erry = max(abs(y_est-yT)./(max(abs(yT))));
%     cond_bad = reshape(not(erry <1),size(y,[2:NLY])); perc_Err = sum(not(erry<1),'all')/numel(erry);
%     if perc_Err > 0.01; fprintf('There are %.2f %% of NaN or infty values\n\n', perc_Err*100); end
%     a(cond_bad) = mean(a,'all','omitnan');  b(cond_bad) = mean(b,'all','omitnan');   c(cond_bad) = mean(c,'all','omitnan');     
%     gaussEqn = 'a*exp(-(1/2)*((x-b)/c)^2)'; startPoints = [mean(a,'all','omitnan'),mean(b,'all','omitnan'),mean(c,'all','omitnan')];
%     logyb = logy(:,cond_bad); Ncond = sum(cond_bad(:)>0); logybm = reshape(logyb.',[Ncond,1,Nx]); xm = reshape(x.',[1,1,Nx]);
%     ad = zeros(Ncond,1); bd = zeros(Ncond,1); cd = zeros(Ncond,1); 
%     for it = 1:Ncond; disp(it); f1 = fit(x,logyb(:,1),gaussEqn,'Start', startPoints); ad(it) = f1.a; bd(it) = f1.b; cd(it) = f1.c; end
%     yb = ad.*exp(-(1./2).*((x'-bd)./cd).^2); yb = yb'; yTb = (yT(:,cond_bad)); err_n(yb,yTb)
%     for it = 1:Ncond; cond_it =  end
        
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


