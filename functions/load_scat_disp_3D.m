% This function calculates the displacement of the scatterers in z-direction after the mechanical loading, in presence of series of strains.
% It is a generalisation of load_scat_Disp.m

% INPUT:
% "alpha": Dimension (Nv,Nalpha_conf). Array of strains. each strain defines a layer. "Nalpha" is hte number of strains.
% "v" Dimension (Nv,3). Array of the "unloaded" z-coordinates
% "varargin{1} = initial_disp". It is the initial bulk displacement of the scatterers. If it is missed then is set as "0", which means that the shallowest border is fixed. This displacement does not depend on the array of strains.

% OUTPUT:
% "v_LT --> dim(Nv,3,Nalpha_conf)" New positon of the scatterers z-coordinates after the meccanical loading of all Nalpha_conf strains in all scatterers
% "varargout{1} = Disp". Dimension(Nalpha+1,1): Array of the borders after the mechanical loading.


function v_LT = load_scat_disp_3D(ALPHA,v,varargin)

if not(size(ALPHA,1) == size(v,1)); error('Strains array must have the same dimension of the scatterers array!'); end

initial_disp = 0;  dc = 5e-6; dcx = dc; dcy = dc; dcz = dc; 
if nargin >= 3; dcx = varargin{1}; end
if nargin >= 4; dcy = varargin{2}; end
if nargin >= 5; dcz = varargin{3}; end
if nargin >= 6; initial_disp = varargin{4}; end

Nalpha_conf = size(ALPHA,2);
x = v(:,1); y = v(:,2); z = v(:,3);
xmin = min(x); xmax = max(x); xbar = xmin:dcx:(xmax+dcx); Nx = length(xbar);
ymin = min(y); ymax = max(y); ybar = ymin:dcy:(ymax+dcy); Ny = length(ybar);
zmin = min(z); zmax = max(z); zbar = zmin:dcz:(zmax+dcz); Nz = length(zbar);
Nxyz = Nx*Ny*Nz;

v_LT = zeros(size(v,1),size(v,2),Nalpha_conf);
for ia = 1:Nalpha_conf
    fprintf('');
    alpha = ALPHA(:,ia);
    v_L = v; z_L = v_L(:,3) + NaN;

    Cx = cell(Nx,1); for ix = 1:(Nx-1); Cx{ix} =( (xbar(ix) <= x) & (xbar(ix+1) > x)); end
    Cy = cell(Ny,1); for iy = 1:(Ny-1); Cy{iy} =( (ybar(iy) <= y) & (ybar(iy+1) > y)); end
    Cz = cell(Nz,1); for iz = 1:(Nz-1); Cz{iz} =( (zbar(iz) <= z) & (zbar(iz+1) > z)); end

    for ix = 1:(Nx-1)
        fprintf('Scatterers displacement estimation: (%d,%d)/(%d,%d)\n',ia,ix,Nalpha_conf,Nx-1)
        for iy = 1:(Ny-1)
            condxy = Cx{ix}&Cy{iy};
            Disp = zeros(Nz,1); Disp(1) = zbar(1) + initial_disp;
            for iz = 1:(Nz-1)
                cond = condxy & Cz{iz};
                sc = sum(cond);
                cond_empty = (sc == 0);
                t = 1.2;
                while cond_empty &  sc < Nxyz;
                    cond = abs(xbar(ix)-x) <= dcx*t & abs(ybar(iy)-y) <= dcy*t & abs(zbar(iz)-z) <= dcz*4*t; sc = sum(cond);
                    t= t*1.5;  cond_empty = (sc == 0);
                end
                if sc == 0; error('There are no strains in the region of scatterers analyzed'); end

                a_mean = mean(alpha(cond)); %alpha_mean(ix,iy,iz) = mean(alpha(cond));
                if isnan(a_mean); a_mean = 0; end
                z_L(cond) = Disp(iz) + (z(cond)-zbar(iz)).*(1-a_mean); % new scatterer positions of the i-th layer after the loading
                Disp(iz+1) = Disp(iz) + dcz.*(1-a_mean); % Displacement of the borders after the loading
            end
        end
    end
    v_L(:,3) = z_L; v_LT(:,:,ia) = v_L;
end
% varargout{1} = Disp;

end
