function [y,varargout] = slit_mean(x,n,varargin);

if nargin == 3
    i0 = varargin{1};
else
    i0 = 1;
end

sx = size(x);
Nx = sx(i0);

N = Nx-n +1; 
n2 = ceil(n/2); 
sy = sx; sy(i0) = N;
y = zeros(sy);

if i0 == 1
    for i = 1:N
        y(i,:,:,:,:,:,:) = sum(x(i:i+n-1,:,:,:,:,:,:),i0)/n;
    end
elseif i0 == 2
   for i = 1:N
        y(:,i,:,:,:,:,:) = sum(x(:,i:i+n-1,:,:,:,:,:),i0)/n;
   end
elseif i0 == 3
   for i = 1:N
        y(:,:,i,:,:,:,:) = sum(x(:,:,i:i+n-1,:,:,:,:),i0)/n;
   end
elseif i0 == 4
   for i = 1:N
        y(:,:,:,i,:,:,:) = sum(x(:,:,:,i:i+n-1,:,:,:),i0)/n;
   end
elseif i0 == 5
   for i = 1:N
        y(:,:,:,:,i,:,:) = sum(x(:,:,:,:,i:i+n-1,:,:),i0)/n;
   end
elseif i0 == 6
   for i = 1:N
        y(:,:,:,:,:,i,:) = sum(x(:,:,:,:,:,i:i+n-1,:))/n;
   end
elseif i0 == 7
   for i = 1:N
        y(:,:,:,:,:,:,i) = sum(x(:,:,:,:,:,:,i:i+n-1))/n;
   end
end


varargout{1} = n2:(n2+N-1);



end

