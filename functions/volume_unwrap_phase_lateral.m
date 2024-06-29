function out_pd_uw = volume_unwrap_phase_lateral(unwrap_filename, phase_matfile, unwrap_config, varargin)
%% out_pd_uw = volume_unwrap_phase_lateral(unwrap_filename, pd_matfile, unwrap_config)
% out_pd_uw = volume_unwrap_phase_lateral(unwrap_filename, pd_matfile, unwrap_config, bscan_extent)
%
% Given a 3-D volume of phase (difference) data, in the form of complex
% numbers, where the amplitude is the weight, and the angle is the
% phase, generates an output volume of unwrapped phase (difference)
%
% INPUT:
% ------
% unwrap_filename: STRING, name of the file to save the unwrapped phase
% information. The unwrapped phase angle data will be saved to the
% dataset 'pd_unwrapped_xyz' in single precision floating point.
%
% phase_matfile: MATFILE, a 'matfile' object used to access the phase
% (difference) data. Expects the complex phase difference to be in the
% dataset 'cplx_phase_diff_xyz'
%
% unwrap_config: STRUCTURE,
%   .PHASE_UNWRAP_Z_RANGE
%   .PHASE_UNWRAP_X_RADIUS
%   .PHASE_UNWRAP_Y_RADIUS
%   .PHASE_UNWRAP_ALWAYS
%   - Axial and lateral ranges used to define the neighbourhood for
%   unwrapping the phase. If UNWRAP_ALWAYS is true, then always unwrap,
%   otherwise unwrap a pixel only if it's weight is less than the
%   weights of the pixels in the neighbourhood.
%
% bscan_extent: ARRAY [start_idx, stop_idx] (OPTIONAL), first to last
% bscan (inclusive) over which the unwrapping should be performed. Used
% if phase_matfile includes less data than its size would otherwise
% indicate.
%
% OUTPUT:
% -------
% out_pd_uw: MATFILE, a 'matfile' object usable for accessing the
% unwrapped phase.
import OCT_OCE.FileIO.*;
import OCT_OCE.OCE.*;

%% Allocate space fo 3-D unwrapped phase difference information
% phase_size = size(phase_matfile, 'cplx_phase_diff_xyz');
oct_max_size_pix.x = size(phase_matfile.cplx_phase_diff_xyz, 1);
oct_max_size_pix.y = size(phase_matfile.cplx_phase_diff_xyz, 2);
oct_max_size_pix.z = size(phase_matfile.cplx_phase_diff_xyz, 3);

oct_vox_size_m = phase_matfile.voxel_size_metres;
out_pd_uw = alloc_file_matfile(unwrap_filename, 'pd_unwrapped_xyz', oct_max_size_pix, oct_vox_size_m, 'single');

%% Determine the bscan ranges
p = inputParser();
p.addOptional('bscan_extent', [1, oct_max_size_pix.y]);
p.parse(varargin{:});

out_first_bscan_idx = p.Results.bscan_extent(1);
out_last_bscan_idx  = p.Results.bscan_extent(2);

% fprintf('Unwrapping %d B-scans: [%d -> %d]\n', ...
%         out_last_bscan_idx - out_first_bscan_idx + 1, ...
%         out_first_bscan_idx, out_last_bscan_idx);

%% Buffer up sufficient en face phase differences to perform unwrapping.
% Initialise the buffer with the first
% unwrap_config.PHASE_UNWRAP_Z_RANGE en face slices. The buffer only
% needs to be big enough to hold unwrap_config.PHASE_UNWRAP_Z_RANGE
% slices.
%
% Only consider the B-scans between the range indicated by bscan_extent
cplx_pd_buffer_xyz = phase_matfile.cplx_phase_diff_xyz(...
    :, out_first_bscan_idx:out_last_bscan_idx, ...
    1:unwrap_config.PHASE_UNWRAP_Z_RANGE);
pd_enface_buffer_xyz = angle(cplx_pd_buffer_xyz);
snr_enface_buffer_xyz = abs( cplx_pd_buffer_xyz);
% The initial en face slices are not unwrapped
out_pd_uw.pd_unwrapped_xyz(:, out_first_bscan_idx:out_last_bscan_idx, ...
                           1:unwrap_config.PHASE_UNWRAP_Z_RANGE) = ...
    pd_enface_buffer_xyz;

% Start unwrapping from the first en face slice after the buffered
% number of slices
rad_pix_uw_xy = [unwrap_config.PHASE_UNWRAP_X_RADIUS, unwrap_config.PHASE_UNWRAP_Y_RADIUS];
for enface_idx = (unwrap_config.PHASE_UNWRAP_Z_RANGE + 1) : oct_max_size_pix.z
%     fprintf('En-face phase difference %d\n', enface_idx);

    % Read in this en face phase difference
    en_face_pd_xy = squeeze( ...
        phase_matfile.cplx_phase_diff_xyz(...
            :, out_first_bscan_idx:out_last_bscan_idx, enface_idx) );

    pd_xy = angle(en_face_pd_xy);
    snr_xy = abs(en_face_pd_xy);

    % Perform the lateral phase unwrapping of this en face slice
    pd_uw_xy = unwrap_phase_diff_lateral(...
        pd_xy, snr_xy, pd_enface_buffer_xyz, snr_enface_buffer_xyz, ...
        rad_pix_uw_xy, unwrap_config.PHASE_UNWRAP_ALWAYS);

    % Save and shuffle the buffers
    pd_enface_buffer_xyz  = cat(3, pd_enface_buffer_xyz(:, :, 2:end),  pd_uw_xy);
    snr_enface_buffer_xyz = cat(3, snr_enface_buffer_xyz(:, :, 2:end), snr_xy);

    % Save the en face slice to disk
    out_pd_uw.pd_unwrapped_xyz(:, out_first_bscan_idx:out_last_bscan_idx, enface_idx) = permute(pd_uw_xy, [1, 2, 3]);
end
