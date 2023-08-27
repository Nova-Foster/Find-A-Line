%Converted from python file by Nova Foster (They/Them) 29/07/2023

NIST_Data = readtable('Lines_eV.csv');

% List of elements to be used in compare, not done in function for efficiency
main_air_elements = {'H', 'He', 'Ar', 'N', 'O'};
all_air_elements = {'H', 'He', 'Ar', 'N', 'O', 'Ne', 'Kr', 'Xe', 'I'};

data = convert_imd("10umVertical - 1000Hz - 950V - 10usmm -ddg - 3.00002 -00028.IMD");
data_subd = background_subtraction(data,700,3);
data_calib = add_calibration(data,700,3,0);

[intensity,wavelengths,times,calib_info] = separate_data_and_calibration(data_calib);
%spark_plot(intensity,wavelengths,times)
%integrated_line_image(intensity,wavelengths)
values = auto_peaks(intensity,wavelengths,true)


function lines = compare(Observed, margin, source)
    % Compare observed wavelength (nm) to NIST values within range +-margin.
    % Source being elements to select from (Main Air, All Air, or All)
    
    % Change format of inputs
    observed_float = str2double(Observed);
    margin_float = str2double(margin);
    source_low = lower(source);
    
    % Select all wavelengths from NIST that fit the range
    lines = NIST_Data(NIST_Data.obs_wl_air_nm >= (observed_float - margin_float) & NIST_Data.obs_wl_air_nm <= (observed_float + margin_float), :);
    
    % Select the sources from that range
    % Main Air: H, He, Ar, N, O
    if strcmp(source_low, 'main air')
        lines = lines(ismember(lines.element, main_air_elements), :);
    % All Air: H, He, Ar, N, O, Ne, Kr, Xe, I
    elseif strcmp(source_low, 'all air')
        lines = lines(ismember(lines.element, all_air_elements), :);
    elseif ~strcmp(source_low, 'all')
        element = {source};
        lines = lines(ismember(lines.element, element), :);
    % Else covers all possible sources
    end
end

function data_2d_array = convert_imd(file_name)
    % Create 2D array from imd: 1392 columns of wavelength & 1040 rows for time
    
    file = fopen(file_name, 'rb'); % Open the file reading as non-text
    
    if file
        % Read first 2 bytes, least important byte first, convert to int from hex
        VER = fread(file, 1, 'uint16', 0, 'l');
        W = fread(file, 1, 'uint16', 0, 'l');
        H = fread(file, 1, 'uint16', 0, 'l');
        
        % Read the rest of the data
        m_s32data = fread(file, W * H, 'int32', 0, 'l');
        
        % Convert the data into floats, with error handling for 0
        fdata = double(m_s32data);
        fdata(m_s32data ~= 0) = fdata(m_s32data ~= 0) / 1000.0;
        
        % Create 2D array
        rows = reshape(fdata, W, H)'; % Reshape data into 2D array
        data_2d_array = flipud(rows);  % Flip data so time 0 is at index 0
        fclose(file);
    else
        data_2d_array = [];
    end
end

function wavelength_whole = cont_wavelengths(calibration)
    % Determine the wavelength value for each pixel based on the calibration file
    % Working based on calibration being 2D array: 1st row = pixel, 2nd row = wavelength
    
    % Find wavelength values per pixel between calibration point 1 and 2
    known_points = size(calibration, 2);
    wavelength_whole = zeros(1, known_points);
    
    current_wavelength = calibration(2, 1);
    next_wavelength = calibration(2, 2);
    pixel_distance = round(calibration(1, 2) - calibration(1, 1));
    
    % Find linearity between first two measurements
    temp_scale = linspace(current_wavelength, next_wavelength, pixel_distance);
    
    % Apply to pixels before first measurement: assumes same linearity as between point 1 and 2
    wavelength_change = temp_scale(2) - temp_scale(1);
    start_val = current_wavelength - (wavelength_change * calibration(1, 1));
    start_scale = linspace(start_val, current_wavelength, round(calibration(1, 1)));
    
    % Create array to be output: currently has wavelength values from pixel 0 to calibration pixel 2
    wavelength_whole = [start_scale, temp_scale];
    
    % Loop for each calibration value
    for i = 1:(known_points - 1)
        try
            current_wavelength = calibration(2, i + 1); % Current wavelength is +1 as the first calibration value has already been used
            next_wavelength = calibration(2, i + 2);
            pixel_distance = round(calibration(1, i + 2) - calibration(1, i + 1)); % Find distance in pixels between the two to be checked
            temp_scale = linspace(current_wavelength, next_wavelength, pixel_distance); % Create the region for current wavelength to next
            
            wavelength_whole = [wavelength_whole, temp_scale]; % Append that to what has already been done
        catch
            break % Exit the loop, should be at the last calibrated value
        end
    end
    
    % Continue the scale for the last value to the final pixel
    wavelength_change = wavelength_whole(end) - wavelength_whole(end - 1); % Amount of wavelength change per pixel
    end_val = next_wavelength + (wavelength_change * (1392 - calibration(1, end))); % Final value of the scale
    end_scale = linspace(calibration(2, end), end_val, 1392 - round(calibration(1, end))); % Create a scale between the final calibration and the end
    wavelength_whole = [wavelength_whole, end_scale]; % Append to pre-existing scale
    
    % Plot the wavelengths vs pixel to show the linearity
    plot(linspace(0, 1391, 1392), wavelength_whole, 'LineWidth', 1);
    xlabel('Pixel');
    ylabel('Wavelength(nm)');
    grid on;
    title('Wavelength Calibration');
    
    % Show the plot
    legend('Wavelength vs Pixel');
    
    % Return the wavelength array
    end

function values = add_calibration(data, center, grating, speed)
    % Add wavelength and speed calibration
    
    % Load correct wavelength file
    wave_filename = sprintf('%d,%d.txt', grating, center);
    raw_wavelengths = load(wave_filename);
    wavelengths = cont_wavelengths(raw_wavelengths);
    calib_info = 'Wave';
    
    % Stack data in the correct way
    values_with_wave = [wavelengths; data];
    
    % Select correct time
    % TODO: Speed is currently just hard coded, need proper timebases in a file and a way to select them
    % all_speeds = load('time_bases.txt');
    % end_time = all_speeds(select based on speed);
    end_time = 0.718 * 200;
    times = linspace(0, end_time, 1040);
    times = [ -1, times ]; % [0,0] of the calibrated file isn't used, -1 to show calibration has been added
    
    values = zeros(1041, 1393);
    for i = 1:1041
        values(i, :) = [times(i), values_with_wave(i, :)]; % Add the time to each column
    end
end

function result = background_subtraction(data, center, grating)
    % Subtract the background
    
    % TODO: Select the background file based on the center and grating
    % back_ground = convert_imd('path_to_background_file.imd');
    back_ground = convert_imd('center 700 Dark.imd');
    
    result = data - back_ground;
end

function [intensity, wavelengths, times, calib_info] = separate_data_and_calibration(data_2d)
    % Separate the whole 2D array into: intensity, wavelength, time, and info character
    
    intensity = data_2d(2:end, 2:end);
    wavelengths = data_2d(1, 2:end);
    times = data_2d(2:end, 1);
    calib_info = data_2d(1, 1);
end

function spark_plot(intensity, wavelengths, times)
    % Plot the image in terms of pixels, horizontal, and vertical profiles
    
    % Color plot of the whole image
    figure;
    imagesc(intensity);
    colormap('hot');
    xlabel('x pixel');
    ylabel('y pixel');
    colorbar;
    
    % Horizontal profile
    figure;
    hori = sum(intensity, 1);
    plot(wavelengths, hori);
    xlabel('Wavelengths');
    ylabel('Intensity');
    
    % Vertical profile
    figure;
    vert = sum(intensity, 2);
    plot(times, vert);
    xlabel('Time (ms)');
    ylabel('Intensity');
end

function integrated_line_image(intensity, wavelengths)
    % Creates a time integrated line image in grayscale
    
    hori = sum(intensity, 1);
    
    % Convert relative data to grayscale values (0 to 1)
    relative_data = hori / max(hori);
    
    figure;
    
    % Create a matrix of size (num_wavelengths, 1000) with relative_data values
    relative_data_matrix = repmat(relative_data, 1000, 1);
    
    % Create an image using imagesc
    imagesc(wavelengths, linspace(0, 1, 1000), relative_data_matrix);
    
    colormap('gray'); % Set grayscale colormap
    
    xlabel('Wavelength (nm)');
    title('Integrated Line Image');
end

function wavelength_over_time(intensity, wavelengths, times, wavelength, index)
    % Plot a wavelength over time
    % Rounds to 3 decimal places
    
    if index == 0 && wavelength ~= 0
        % Select the wavelength based on calibration
        index = find(wavelengths == wavelength);
    end
    
    % Select the correct intensity
    specific_intensity = intensity(:, index);
    
    % Plot and set the correct title
    over_time = figure;
    plot(times, specific_intensity);
    xlabel('Time (ms)');
    ylabel('Intensity');
    title_str = sprintf('%.3f over time', round(wavelength, 3));
    title(title_str);
    
    hold off;
end

function [pks_intensity,pks_xs] = auto_peaks(intensity,wavelengths,strict)
    figure;
    hori = sum(intensity, 1);
    plot(wavelengths, hori);
    xlabel('Wavelengths');
    ylabel('Intensity');

    if strict ==false
    rel_prom=single(0.075);
    rel_width=single(0.2);
    else
    rel_prom = single(0.15);
    rel_width=single(0.25);
    end
    
  
    pks = findpeaks(hori,"Annotate","peaks","MinPeakProminence",rel_prom,"MinPeakWidth",rel_width);
    findpeaks(hori,"MinPeakProminence",rel_prom,"MinPeakWidth",rel_width)

end