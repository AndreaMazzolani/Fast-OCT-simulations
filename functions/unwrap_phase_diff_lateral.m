function pd_uw_xy = unwrap_phase_diff_lateral(pd_xy, snr_xy, prev_pd_xyz, prev_snr_xyz, rad_pix_uw_xy, UNWRAP_ALWAYS)
%% pd_uw_xy = unwrap_phase_diff_lateral(pd_xy, snr_xy, prev_pd_xyz, prev_snr_xyz, rad_pix_uw_xy, UNWRAP_ALWAYS)
% Performs 3-D lateral phase unwrapping for optical coherence
% elastography. Unwraps the phase difference in the en-face plane given
% the phase differences in the n en-face planes preceeding it.
%
% INPUT:
% ------
%
% pd_zx: (x, y) ordering. This is the input phase difference to unwrap,
% and will be unwrapped based on the previous (in z) phase differences
% given by prev_pd_zxy. These previous (in z) phase differences are all
% assumed to be already unwrapped.
%
% snr_xy: (x, y) ordering. Associated OCT signal-to-noise ratio of
% each phase difference measurement. This should be in linear scale.
%
% prev_pd_xyz: (x, y, z) ordering. The previously (in z) unwrapped phase
% differences used to unwrap the current en face plane.
%
% prev_snr_xyz: (x, y, z) ordering. The OCT SNR associated with each of
% the prev_pd_zxy measurements.
%
% rad_pix_uw_xy: [x-pix, y-pix]. Radius of pixels in the x- and
% y-dimensions to use when unwrapping. This will look in a neighbour of
% rad_pix_uw_xy before and after the current pixel.
%
% UNWRAP_ALWAYS: true or false. If true, always unwrap the current pixel
% based on its neighbourhood. If false, then the current pixel is only
% unwrapped if the weights of its neighbourhood out weigh the weight
% associated with the current pixel.
%
% OUTPUT:
% -------
%
% pd_uw_xy: (x, y) ordering. Phase difference B-scan unwrapped.

%% Check the inputs
if numel(rad_pix_uw_xy) ~= 2
    error('%s: rad_pix_uw_xy: expected 2 elements!', mfilename() );
end

% Unwrapping the en face slice based on the previous en face slices
tmp_enface_x_y     = pd_xy;
weights_enface_x_y = snr_xy;

%% Axial Unwrapping BEGIN
prev_pd_x_y_z     = prev_pd_xyz;
prev_weight_x_y_z = prev_snr_xyz;

%% Variance of a weighted mean =
% var_mean = 1 / (sum(w_i)),
% where w_i = 1 / var_i, var_i = variance of the i-th value
% => weight_mean = sum(w_i)
avg_pd_weight_x_y = sum( prev_weight_x_y_z, 3, 'double' );
% Sanity check - since the weight is a sum of SNRs, anything less than
% the number of elements (min SNR = 1, so min sum of n SNRs = n) is
% unreliable. In particular, we can set the minimum sum of SNRs to 1, as
% the only time its less than this will be terribleness in the noise, or
% a weight of 0 due to incompletely processed data.
% avg_pd_weight_x_y(avg_pd_weight_x_y < 1) = 1;

%% Calculate the weighted mean phase difference at each (x, y)
%% location, of the previous unwrap_z_pix en-face scans
avg_pd_x_y = sum(prev_pd_x_y_z .* prev_weight_x_y_z, 3, 'double') ...
    ./ avg_pd_weight_x_y;

%% Unwrap every pixel in the current en-face scan based on the
%% averaged previous PD of the preceeding unwrap_z_pix scans
% look at the difference between the averaged previous PD, and the
% current PD
diff_pd_x_y = tmp_enface_x_y - avg_pd_x_y;

if (UNWRAP_ALWAYS)
    % Unconditionally unwrap the current pixel based on the
    % self.unwrap_z_pix pixels above it.
    closest_x_y = round( diff_pd_x_y / (2*pi) );
else
    % Initialise the unwrapping array to 0, meaning by default, don't unwrap
    closest_x_y = zeros(size(tmp_enface_x_y), class(tmp_enface_x_y));

    % But only if the current pixel is less trustworthy than the previous
    % window (ie, the weight of the current pixel < weight of the
    % previous window), then calculate how many multiples of 2*pi do
    % we need to minimise this difference? Or, what's the closest
    % multiple of 2*pi to this difference?
    closest_x_y( weights_enface_x_y <= avg_pd_weight_x_y ) = ...
        round( diff_pd_x_y( weights_enface_x_y <= avg_pd_weight_x_y ) / (2*pi) );
end

% Subtract this closest multiple of 2*pi from the current phase
% difference to yield the 'unwrapped' phase difference
tmp_enface_x_y = tmp_enface_x_y - closest_x_y * 2*pi;

%% Axial Unwrapping COMPLETE

%% Lateral Unwrapping BEGIN
% Only laterally unwrap if the unwrapping radius is large enough that it makes sense
if any(rad_pix_uw_xy) > 0
    % At every (x, y) location on the current En-Face slice, unwrap the
    % pixel based on the neighbourhood of pixels within a lateral window
    % within rad_pix_uw_xy pixels distance.
    %
    % First generate the convolution kernel representing the
    % neighbourhood around each (x, y) pixel.
    %
    % eg. if rad_pix_uw_xy = [1, 1]:
    % kern = [ 1 1 1
    %          1 0 1
    %          1 1 1 ];
    neigh_sum_kern_xy = ones( 2 * rad_pix_uw_xy + 1 );
    neigh_sum_kern_xy(rad_pix_uw_xy(1) + 1, rad_pix_uw_xy(2) + 1) = 0;

    % Calculate the weights of the phase difference pixels in the
    % neighbourhood around each (x, y) location
    neigh_weights_x_y = conv2(double(weights_enface_x_y), ...
                              double(neigh_sum_kern_xy), 'same');
    % Sanity check - as above, we can set the minimum sum of SNRs to 1, as
    % the only time its less than this will be terribleness in the
    % noise, or a weight of 0 due to incompletely processed data.
    % neigh_weights_x_y(neigh_weights_x_y < 1) = 1;

    % Calculate the weighted average phase difference of the pixels in the
    % neighbouhood around each (x, y) location. This is based off the
    % axially unwrapped phase-difference.
    neigh_avg_pd_x_y = conv2(...
        double(tmp_enface_x_y) .* double(weights_enface_x_y), ...
        double(neigh_sum_kern_xy), 'same') ...
        ./ neigh_weights_x_y;

    % Get the difference between the phase difference of the neighbourhood,
    % and the phase difference of the pixel in each (x, y)
    neigh_diff_pd_x_y = tmp_enface_x_y - neigh_avg_pd_x_y;

    if (UNWRAP_ALWAYS)
        % Unconditionally unwrap the current pixel based on the neighbourhood
        % around it.
        neight_closest_x_y = round( neigh_diff_pd_x_y / (2*pi) );
    else
        % Initialise the unwrapping array to 0, meaning by default, don't unwrap
        neight_closest_x_y = zeros(size(tmp_enface_x_y), class(tmp_enface_x_y));

        % For every (x, y) pixel, if the weight of the pixel is less than the
        % weights of the combined neighbouring pixels, then unwrap the pixel
        % at (x, y) by adding a multiple to 2*pi to minimise the difference
        % between the pixel at (x, y) and the weighted mean of the pixels in
        % the neighbourhood around (x, y)
        neight_closest_x_y( weights_enface_x_y <= neigh_weights_x_y ) = ...
            round( neigh_diff_pd_x_y( weights_enface_x_y <= neigh_weights_x_y ) / (2*pi) );
    end

    % Subtract this closest multiple of 2*pi from the current phase
    % difference to yield the 'LATERALLY unwrapped' phase difference
    tmp_enface_x_y = tmp_enface_x_y - neight_closest_x_y * 2*pi;
end
%% Lateral Unwrapping COMPLETE

pd_uw_xy = tmp_enface_x_y;
