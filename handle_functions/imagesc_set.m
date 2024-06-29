function varargout = imagesc_set(B, varargin)

sB = size(B);

if length(size(B)) > 2; B = squeeze(B); if length(size(B)) > 2; error('The dataset dimension is too high!'); end; end;

sB = size(B);
Nz = sB(1); 
Nx = sB(2);

nofig = 0;
idx = 0;

if nargin == 1 
    z_vec = 1:Nz; 
    x_vec = 1:Nx;
    ns = 0.05;
elseif nargin == 2
    z_vec = varargin{1}; 
    x_vec = 1:Nx;
    ns = 0.05;
elseif nargin == 3
    z_vec = varargin{1};
    x_vec = varargin{2};
    ns = 0.05;
elseif nargin == 4
    z_vec = varargin{1};
    x_vec = varargin{2};
    idx = varargin{3};
elseif nargin == 5
    z_vec = varargin{1};
    x_vec = varargin{2};
    idx = varargin{3};
    ns = varargin{4};
else 
    error('Too many input arguments!');
end

if not(length(z_vec) == Nz) 
    if (length(z_vec) == Nx) & (length(x_vec) == Nz)
        B = B.';
    else
        error('Dimensions of the arrays do not match!');
    end
end

% if not( (length(z_vec) == Nz) & (length(x_vec) == Nx) & (length(ns) == 1) & (length(idx) == 1) ) 
%     error('Dimensions of the arrays do not match!');
% end

if not( (length(ns) == 1) & (length(idx) == 1) ) 
    error('Dimensions of the arrays do not match!');
end

[m_low,m_sup] = fig_set(B, ns);

if not(idx == 0); figure(idx); clf; end
imagesc(x_vec,z_vec,B); caxis([m_low,m_sup]); colorbar;
varargout{1} = m_low;
varargout{2} = m_sup;


end