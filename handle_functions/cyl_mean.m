function [Mat_mean,varargout] = cyl_mean(Mat,Ir_c,Ix_c,varargin)

% This function calculates the mean of the multi-dimensionale array "Mat"
% in the indices related to a 2D sphere or 3D cylinders with radius Ir_c
% and center [Ix_c,Iy_c] with z indices "Indz"

sM = size(Mat); Lm = length(sM); if Lm == 2 & min(sM) == 1; Lm = 1; end
if not(Lm == nargin-2); error('The dimensions of the array and the box does not match'); end

IX = reshape(1:sM(1),[sM(1),1,1]); IY = 1; Iy_c = 1; Indz = 1;

if nargin >= 4; Iy_c = varargin{1}; IY = reshape(1:sM(2),[1,sM(2),1]); 
    if Lm == 3; Indz = reshape(1:sM(3),[1,1,sM(3)]); end;  % nargin == 4 the cylinder is taken entirely
end
if nargin >= 5; Indz = varargin{2};   end
if nargin >= 6; error('Arrys having dimensions higher than 5 are not included'); end

Mat_round_condxy = ((IX-Ix_c).^2 + (IY-Iy_c).^2 < Ir_c.^2);
Mat_round_condz  = Mat.*0; Mat_round_condz(:,:,Indz) = 1;
Mat_round_cond = Mat_round_condz.*Mat_round_condxy;

clear Mat_round_condxy Mat_round_condz
Mat_mean = sum(Mat_round_cond.*Mat,'all')./sum(Mat_round_cond(:)); 

if nargout == 2; varargout{1} = Mat_round_cond; end

end