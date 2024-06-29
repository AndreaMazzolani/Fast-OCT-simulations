function [Mat_mean,varargout] = box_mean(Mat,Indx,varargin)

% This function calculates the mean of the multi-dimensionale array "Mat" in the indices box --> [Indx,Indy,Indz,Indt,Indw]

Indy = 1; Indz = 1; Indt = 1; Indw = 1;
if nargin >= 3; Indy = varargin{1}; end
if nargin >= 4; Indz = varargin{2}; end
if nargin >= 5; Indt = varargin{3}; end
if nargin >= 6; Indw = varargin{4}; end
if nargin >= 7; error('Arrys having dimensions higher than 5 are not included'); end

sM = size(Mat); Lm = length(sM); if Lm == 2 & min(sM) == 1; Lm = 1; end
if not(Lm == nargin-1); error('The dimensions of the array and the box does not match'); end


Mat_box_cond = Mat.*0; Mat_box_cond(Indx,Indy,Indz,Indt,Indw) = 1;
Mat_mean = sum(Mat_box_cond.*Mat,'all')./sum(Mat_box_cond(:));

if nargout == 2; varargout{1} = Mat_box_cond; end

end