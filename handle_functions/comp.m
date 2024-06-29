
function comp(y0,varargin);

sy0 = size(y0);
Ly = length(sy0);

if Ly > 2 | not(max(sy0) == prod(sy0))
    error('The functions to compare must be arranged in a 1D or 2D array!');
end

Y = y0(:);
condx = 0;

for in = 1:(nargin-1)
    if isempty(varargin{in});
        condx = 1; ix = 0; num = 1;
    else
        if condx == 0
            Y(:,in+1) = varargin{in}(:);
        else
            if ix == 0;
                x = varargin{in};
                ix = ix+1;
            elseif ix == 1;
                num = varargin{in};
            end
        end
    end
end


sy = size(Y);
if condx == 0;
    x = 1:sy(1); num = 1;
else
    if ix == 0;
        num = 1;
    end
end

Lx = length(x);

if not(Lx == sy(1)); error('The size of domain and codmain must match!'); end

err_m  = @(f1,f2) 2.*abs(f1(:) - f2(:))./(abs(f1(:))+abs(f2(:)));

figure(num); clf;
if max(abs(imag(Y(:)))) == 0;
    if min(sy) == 1;
      y = Y;  plot(x,y);
    else;
        for k = 1:sy(2)
            y = Y(:,k);
            subplot(2,1,1); hold on;  plot(x,y);  title('funcs');
            if k > 1
                subplot(2,1,2); hold on;  plot(x,err_m(y(:),y_old(:)));  title('rel. errs.');
            end
            y_old = y;
        end
    end

else
    if min(sy) == 1;
        y = Y;
        subplot(3,1,1); hold on; plot(x,real(y) );  title('real');
        subplot(3,1,2); hold on; plot(x,imag(y) ); hold on; title('imag');
        subplot(3,1,3); hold on; plot(x,abs(y) );  hold on; title(' abs');
    else
        for k = 1:sy(2)
            y = Y(:,k);
            %     if nargin > 3 ; laby = [num2str(k)]; else; laby = varargin{k}(); end;
            subplot(3,3,1); hold on; plot(x,real(y));
            subplot(3,3,4); hold on; plot(x,imag(y));
            subplot(3,3,7); hold on; plot(x, abs(y));

            if k > 1
                subplot(3,3,2); hold on; plot(x,err_m(real(y),real(y_old)));
                subplot(3,3,5); hold on; plot(x,err_m(imag(y),imag(y_old)));
                subplot(3,3,8); hold on; plot(x, err_m(abs(y),abs(y_old)));
            end
            subplot(sy(2),3,3*k); hold on; plot(x,abs(y));

            y_old = y;
        end

        subplot(3,3,1); hold on;  ylabel('Real'); title('Arrays');
        subplot(3,3,4); hold on; ylabel('Imag'); title('imag');
        subplot(3,3,7); hold on; ylabel('Abs'); title(' abs');
        subplot(3,3,2); hold on; title('Subsequent relative errors');
        subplot(3,3,5); hold on;
        subplot(3,3,8); hold on;
        subplot(sy(2),3,3); hold on; title(' Abs of single  Arrays');
        sgtitle('Arrays comparison');
    end
end

end



