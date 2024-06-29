function [strain_opt,varargout] = strain_optimized(strain,error,Nx_slit,Ny_slit,Nz_slit,varargin)
% The strain and error matrices must have coordinate indices ordered as "z,x,y"


sx = size(strain); IZ = 1:sx(1); IX = 1:sx(2); IY = 1;
if length(sx) > 2; IY = 1:sx(3); end 

Nerr = 12;
if nargin >= 6; Nerr = varargin{1}; end
if nargin >= 7; error('Tthere are too many input data' ); end

% str_tmp = squeeze(strain); reg2 = real(1./( squeeze(error.^Nerr)));
str_tmp = strain; reg2 = real(1./( error.^Nerr));

if Nx_slit > 1; [reg2_x,IX] = slit_mean(reg2,Nx_slit,2); str_optx = slit_mean(str_tmp.*reg2,Nx_slit,2)./reg2_x;
else; reg2_x = reg2;str_optx = str_tmp; end;

if Ny_slit > 1; [reg2_xy,IY] = slit_mean(reg2_x,Ny_slit,3); str_optxy = slit_mean(str_optx.*reg2_x,Ny_slit,3)./reg2_xy;
else; reg2_xy = reg2_x; str_optxy = str_optx; end;

if Nz_slit > 1; [reg2_xyz,IZ] = slit_mean(reg2_xy,Nz_slit);   strain_opt =  slit_mean(str_optxy.*reg2_xy,Nz_slit)./reg2_xyz;
else; strain_opt = str_optxy; end;

varargout{1} = IZ; varargout{2} = IX; varargout{3} = IY;
 

end





