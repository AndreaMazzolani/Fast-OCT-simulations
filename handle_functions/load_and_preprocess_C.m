function [spectral_data_kx, k_vec, spectrum] = load_and_preprocess_C(filename,index, flatten, apod_spectrum_from_bscan_index)
    % first parameter can be file handle, if the file was opened before


     
    if(isstruct(filename))
        dont_close = true;
        oct_file = filename;
    else
        dont_close = false;
        oct_file = OCTFileOpen(filename);  
        ascan_avg = str2double(oct_file.head.Acquisition.IntensityAveraging.AScans.Text);
        Nr = size(raw_tmp,2)/ascan_avg;
    end

    if( ~exist('flatten','var') || isempty(flatten) )
        flatten=true;
    end
    if( ~exist('apod_spectrum_from_bscan_index','var') || isempty(apod_spectrum_from_bscan_index) )
        apod_spectrum_from_bscan_index=false;
    end
    NB = OCTFileGetNrRawData(oct_file);
    
    if(isnumeric(index))
        if(Nx <= index)
           error("Wrong B scan index");
        end
        [raw,spectrum] = OCTFileGetRawData(oct_file,index);
        
    elseif(index == "avg")
        [raw,spectrum] = OCTFileGetRawData(oct_file,0);

        for ind=1:(NB-1)
            [raw_tmp,spectrum_tmp] = OCTFileGetRawData(oct_file,ind);
            raw = raw + raw_tmp;
            spectrum = spectrum + spectrum_tmp;
        end

        raw = raw / NB;
        spectrum = spectrum / NB;
    elseif(index == "all")
        [raw_tmp,spectrum] = OCTFileGetRawData(oct_file,0);
        raw_kxy = zeros(size(raw_tmp,1),size(raw_tmp,2),NB);
        raw_kxy(:,:,1) = raw_tmp;

    
        for ind=1:(NB-1)
            [raw_tmp,spectrum_tmp] = OCTFileGetRawData(oct_file,ind);
            raw_tmp = reshape(raw_tmp, [numel(spectrum), ascan_avg, Nr  ]);
            raw_tmp = squeeze(mean(raw_tmp, 2));
            raw_kxy(:,:,ind+1) = raw_tmp;
            spectrum = spectrum + spectrum_tmp;
        end

        spectrum = spectrum / NB;
    else
        error("Wrong B scan index");
    end
    
    % average ascan_avg consecutive A-scans into one
    if ascan_avg >= 2
        raw = reshape(raw, [numel(spectrum), ascan_avg, size(raw,2)/ascan_avg  ]);
        raw = squeeze(mean(raw, 2));
    end
    
    % override spectrum with custom spectrum
    % useful if reference arm was closed
    if isnumeric(apod_spectrum_from_bscan_index)
        if isscalar(apod_spectrum_from_bscan_index)
            if(OCTFileGetNrRawData(oct_file) <= apod_spectrum_from_bscan_index)
               error("Wrong apod spectrum B scan index");
            end
            [~,spectrum] = OCTFileGetRawData(oct_file,apod_spectrum_from_bscan_index);
        else
           % use apod_spectrum_from_bscan_index as spectrum
           if size(raw,1) ~= numel(apod_spectrum_from_bscan_index)
                error("Custom apod spectrum has wrong dimensions");
           end
           spectrum = apod_spectrum_from_bscan_index(:);
        end
    end

    
    if flatten
        flat = raw ./ spectrum;
        flat = flat - mean(flat,1);
    else
        flat = raw;
    end
    
    chirp_data = OCTFileGetChirp(oct_file);
%     interpolatingFunction = @(asc) interp1(chirp_data,asc,0:2047,'pchip');
%     spectral_data_kx = interpolatingFunction(flat);
%     spectrum = interpolatingFunction(spectrum);
    spectral_data_kx = flat; %interpolatingFunction(flat);
    
    % flip order as data is decreesing in k originally
    % (see DETERMINATION OF THE WAVELENGTH AXIS FOR A SPECIFIED OCT SYSTEM)
    spectral_data_kx = flip(spectral_data_kx,1);
    spectrum = flip(spectrum);
    
    if(~dont_close)
        OCTFileClose(oct_file);
    end
    
    %the values for the UCL system are:
    %Lambda min: 1170.5 nm
    %Lambda max: 1407.8 nm
    lam_min = 1170.5e-9;
    lam_max = 1407.8e-9;
    
    k_min = 2*pi/lam_max;
    k_max = 2*pi/lam_min;

    k_vec = linspace(k_min,k_max,2048);
end