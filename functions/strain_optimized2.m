function strain_opt = strain_optimized2(strain,error,Nx_slit,Ny_slit,Nz_slit,varargin)
% The strain and error matrices must have coordinate indices ordered as "z,x,y"

sx = size(strain);% IZ = 1:sx(1); IX = 1:sx(2); IY = 1:sx(3);

Nerr = 12; str_nan = '';
if nargin >= 6; Nerr = varargin{1}; end
if nargin >= 7; str_nan = varargin{2}; end
if nargin >= 8; error('There are too many input data' ); end

str_tmp = strain; reg2 = real(1./( error.^Nerr));

if not(strcmp(str_nan,'omitnan'));
    if Nx_slit > 1; reg2_x = movmean(reg2,Nx_slit,2); str_optx = movmean(str_tmp.*reg2,Nx_slit,2)./reg2_x; else; reg2_x = reg2; str_optx = str_tmp; end;
    if Ny_slit > 1; reg2_xy = movmean(reg2_x,Ny_slit,3); str_optxy = movmean(str_optx.*reg2_x,Ny_slit,3)./reg2_xy; else; reg2_xy = reg2_x; str_optxy = str_optx; end;
    if Nz_slit > 1; reg2_xyz = movmean(reg2_xy,Nz_slit);   strain_opt =  movmean(str_optxy.*reg2_xy,Nz_slit)./reg2_xyz; else; strain_opt = str_optxy; end;

else;
    if Nx_slit > 1; reg2_x = movmean(reg2,Nx_slit,2,'omitnan'); str_optx = movmean(str_tmp.*reg2,Nx_slit,2,'omitnan')./reg2_x; else; reg2_x = reg2; str_optx = str_tmp; end;
    if Ny_slit > 1; reg2_xy = movmean(reg2_x,Ny_slit,3,'omitnan'); str_optxy = movmean(str_optx.*reg2_x,Ny_slit,3,'omitnan')./reg2_xy; else; reg2_xy = reg2_x; str_optxy = str_optx; end;
    if Nz_slit > 1; reg2_xyz = movmean(reg2_xy,Nz_slit,'omitnan');   strain_opt =  movmean(str_optxy.*reg2_xy,Nz_slit,'omitnan')./reg2_xyz; else; strain_opt = str_optxy; end;

end

% fprintf('There are %d NAN on %d data\n', sum(isnan(strain_opt(:))),numel(strain_opt));
end





