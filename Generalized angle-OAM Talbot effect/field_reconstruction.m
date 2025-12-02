% A script to reconstruct a complex field from an intefrerogram picture
% Last updated 10.10.2024
% Matias Eriksson (matias.eriksson@tuni.fi)

%% Constants

% For finding the center of mass of filtered Fourier space and real space
find_center_of_masses_flag = 0;
center_of_mass_threshold = 0.2;

% Whether to plot the OAM spectral phase
plot_phases = 0;
filter_low_intensity_phases = 1; % Whether to not display spectral phase for low intensity modes

% Whether to rotate the real field (rotation angle defined later)
rotate_field_flag = 1;

% Whether to correct the linear phase in OAM spectrum due to rotation
correct_linear_phase = 0;

fft_plot_area = 20; % 1/2 of side length of square area of fft to plot

%% Load data & specify settings

% SETTING EXPLANATIONS FOR DATASETS:

% img: off-axis interferogram
% signal_img: signal without reference beam
% cx_fft, cy_fft: shift center of Fourier space to these coordinates
% cx_f, cy_f: shift center of real space to these coordinates
% fft_filter_radius: radius of circular aperture for Fourier filtering
% radi, rado: inner and outer radius of ring for real space filtering
%             (Better OAM spectrum reconstruction, reduce effect of stray  
%             light not in ring-core
% padding_flag: Whether to pad the image with zeroes (effect: interpolation
%               in Fourier space)
% padding_factor: Factor of increasing image size with padding
% cx_fft_offset, cy_fft_offset: Offset the cx_fft,cy_fft in interpolated 
%                               Fourier space
% rotation_angle: Rotation angle of the field in degrees



% UNCOMMENT THE SELECTED DATASET SETTINGS FOR EVALUATION

% Settings 1, 3-to-4
% img = double(imread('.\3-to-4\stage3.png')); % Without 2nd Talbot phase mask
img = double(imread('.\3-to-4\stage4.png')); % With 2nd Talbot phase mask
signal_img = double(imread('.\3-to-4\signal.png'));
cx_fft = 720;
cy_fft = 530;
cx_f = 570;
cy_f = 455;
fft_filter_radius = 0.0125;
radi = 0.25;
rado = 0.5;
padding_flag = 1;
padding_factor = 7;
cx_fft_offset = -1;
cy_fft_offset = 1;
rotation_angle = 86.6;


% % Settings 2, 4-to-3
% img = double(imread('.\4-to-3\stage3.png')); % Without 2nd Talbot phase mask
% % img = double(imread('.\4-to-3\stage4.png')); % With 2nd Talbot phase mask
% signal_img = double(imread('.\4-to-3\signal.png'));
% cx_fft = 718;
% cy_fft = 529;
% cx_f = 570;
% cy_f = 455;
% fft_filter_radius = 0.0125;
% radi = 0.25;
% rado = 1;
% padding_flag = 1;
% padding_factor = 10;
% cx_fft_offset = 1;
% cy_fft_offset = 0;
% rotation_angle = -123.1;


% % Settings 4, 4-to-6 
% img = double(imread('.\4-to-6\stage3.png'));  % Without 2nd Talbot phase mask
% % img = double(imread('.\4-to-6\stage4.png'));  % With 2nd Talbot phase mask
% signal_img = double(imread('.\4-to-6\signal.png'));
% cx_fft = 718;
% cy_fft = 529;
% cx_f = 572;
% cy_f = 455;
% fft_filter_radius = 0.0125;
% radi = 0.25;
% rado = 0.5;
% padding_flag = 1;
% padding_factor = 9;
% cx_fft_offset = 1;
% cy_fft_offset = 1;
% rotation_angle = -123.5;



% % Settings 5, 6-to-4
% % img = double(imread('.\6-to-4\stage3.png')); % Without 2nd Talbot phase mask
% img = double(imread('.\6-to-4\stage4.png')); % With 2nd Talbot phase mask
% signal_img = double(imread('.\6-to-4\signal.png'));
% cx_fft = 718;
% cy_fft = 529;
% cx_f = 570;
% cy_f = 453;
% fft_filter_radius = 0.0125;
% radi = 0.25;
% rado = .5;
% padding_flag = 1;
% % padding_factor = 7;
% % cx_fft_offset = 1;
% % cy_fft_offset = 0;
% % padding_factor = 10;
% % cx_fft_offset = 2;
% % cy_fft_offset = 0;
% padding_factor = 12;
% cx_fft_offset = 3;
% cy_fft_offset = 0;
% rotation_angle = -7.1;



% % Settings 6, 2-to-6
% % img = double(imread('.\2-to-6\stage3.png')); % Without 2nd Talbot phase mask
% img = double(imread('.\2-to-6\stage4.png')); % With 2nd Talbot phase mask
% signal_img = double(imread('.\2-to-6\signal.png'));
% cx_fft = 749;
% cy_fft = 524;
% cx_f = 520;
% cy_f = 453;
% fft_filter_radius = 0.0125;
% radi = 0.25;
% rado = 0.5;
% padding_flag = 1;
% padding_factor = 10;
% cx_fft_offset = 1;
% cy_fft_offset = 1;
% add_global_phase = 3.0;
% rotation_angle = -64.5;


% % Settings 7, 6-to-2
% img = double(imread('.\6-to-2\stage3_4.png')); % The petals have the same phase already after RCF
% signal_img = double(imread('.\6-to-2\signal.png'));
% cx_fft = 749;
% cy_fft = 523;
% cx_f = 520;
% cy_f = 453;
% fft_filter_radius = 0.0125;
% radi = 0.25;
% rado = 0.5;
% padding_flag = 1;
% padding_factor = 6;
% cx_fft_offset = 1;
% cy_fft_offset = 2;
% % padding_factor = 12;
% % cx_fft_offset = 2;
% % cy_fft_offset = 4;
% % padding_factor = 10;
% % cx_fft_offset = 4;
% % cy_fft_offset = -9;
% % cx_fft_offset = 2;
% % cy_fft_offset = 0;
% correct_linear_phase = 0;
% rotation_angle = -72.3;



% % Settings 8, 12-to-2
% % img = double(imread('.\12-to-2\stage3.png')); % Without 2nd Talbot phase mask
% img = double(imread('.\12-to-2\stage4.png')); % With 2nd Talbot phase mask
% signal_img = double(imread('.\12-to-2\signal.png'));
% cx_fft = 749;
% cy_fft = 521;
% % cx_fft = 748;
% % cy_fft = 522;
% cx_f = 518;
% cy_f = 451;
% fft_filter_radius = 0.0125;
% padding_flag = 1;
% % padding_factor = 5;
% % cx_fft_offset = 7;
% % cy_fft_offset = 3;
% % cx_fft_offset = 2;
% % cy_fft_offset = 8;
% % cx_fft_offset = 7;
% % cy_fft_offset = 10;
% padding_factor = 10;
% cx_fft_offset = 10;
% cy_fft_offset = 4;
% % cx_fft_offset = 13;
% % cy_fft_offset = 5;
% radi = 0.25;
% rado = 0.5;
% rotation_angle = -61.2;


% % % Settings 9, 12-to-3
% img = double(imread('.\12-to-3\stage3.png')); % Without 2nd Talbot phase mask
% % img = double(imread('.\12-to-3\stage4.png')); % With 2nd Talbot phase mask
% signal_img = double(imread('.\12-to-3\signal.png'));
% cx_fft = 749;
% cy_fft = 524;
% cx_f = 505;
% cy_f = 455;
% fft_filter_radius = 0.0125;
% padding_flag = 0;
% padding_factor = 10;
% cx_fft_offset = 0;
% cy_fft_offset = 0;
% radi = 0.25;
% rado = 1;
% correct_linear_phase = 1;
% rotation_angle = -182.6;


% % Settings Input, before RCF. Choose any interferogram from the input folder!
% img = double(imread('.\input\OAM modes\l10.png'));
% % img = double(imread('.\input\Petals in phase\gen_12_inphase_400.png'));
% % img = double(imread('.\input\Petals with talbot(or custom) phase\gen_4_(-3_12)_200.png'));
% % signal_img = double(imread('.\input\signal (intensity only)\intensity_ring_445_30.png'));
% 
% cx_fft = 743;
% cy_fft = 536;
% cx_f = 1750;
% cy_f = 1447;
% fft_filter_radius = 0.015;
% padding_flag = 1;
% padding_factor = 10;
% cx_fft_offset = 1;
% cy_fft_offset = 0;
% radi = 0.25;
% rado = 0.5;

%% Initialize

[ny,nx] = size(img); 

% Pad image with zeroes (effect: interpolate Fourier space)
if padding_flag
    % New image size
    ny = ny*padding_factor;
    nx = nx*padding_factor;
    
    % Pad image with zeroes
    img = padortrim_array(img,ny,nx);
    
    % New coordinates to shift the Fourier space center to
    cx_fft = cx_fft*padding_factor+cx_fft_offset;
    cy_fft = cy_fft*padding_factor+cy_fft_offset;
end

% Grids
x=((-nx)/2:(nx)/2-1)/nx;            
y=((-ny)/2:(ny)/2-1)/ny; % Division by nx here would be correct to ensure equal step size in x and y. 
                         % However for Fourier filtering this is better, as Fourier space coordinates have a 
                         % different step size in kx and ky due to non-square image. This division is corrected later in
                         % 'Crop back to original size' section
[X,Y]=meshgrid(x,y);
Rad=sqrt(X.^2+Y.^2);
Angle=angle(X+1i.*Y)+pi; %Matrix with all Angles starting left-center


%% Field reconstruction

% FFT of image
fft_img = fft2(img);

% Filter out zeroth component (easier to visualize fft)
% fft(1,1) = 0;

% Shift zeroth component to middle
fft_img = fftshift(fft_img);

% shift FFT (center 1st diffraction order)
fft_shift = imtranslate(fft_img,[nx/2-cx_fft,ny/2-cy_fft]);

% filter FFT (filter out other diffraction orders)
mask = Rad<fft_filter_radius;
fft_shift_filt = fft_shift.*mask;

% IFFT (field reconstruction)
field = ifft2(ifftshift(fft_shift_filt));

% Visualize filtered and shifted fft
figure(10)
if padding_flag
    imagesc(abs(fft_shift_filt(ny/2-fft_plot_area*padding_factor:ny/2+fft_plot_area*padding_factor,nx/2-fft_plot_area*padding_factor:nx/2+fft_plot_area*padding_factor)));
    hold on
    scatter(fft_plot_area*padding_factor+1,fft_plot_area*padding_factor+1, 'rx');
    hold off
else
    imagesc(abs(fft_shift_filt(ny/2-fft_plot_area:ny/2+fft_plot_area,nx/2-fft_plot_area:nx/2+fft_plot_area)));
        hold on
    scatter(fft_plot_area+1,fft_plot_area+1, 'rx');
    hold off
end
axis square
title('Shifted and filtered Fourier space')


% % Visualize reconstructed field
% figure(11)
% intensity_phase_plot(field);


%% Find center of mass of cropped FFT (for checking)
if find_center_of_masses_flag
    [cx,cy] = center_of_mass(fft_shift_filt,center_of_mass_threshold);
    d_cx = cx-nx/2 % Difference to currently chosen center coordinate
    d_cy = cy-ny/2 % Difference to currently chosen center coordinate 
end

%% Crop back to original size

if padding_flag
    nx_pad = nx;
    ny_pad = ny;
    nx = nx/padding_factor;
    ny = ny/padding_factor;
    
    nx_diff = nx_pad-nx;
    ny_diff = ny_pad-ny;
    
    field = field(1+ny_diff/2:end-ny_diff/2,1+nx_diff/2:end-nx_diff/2);
end

x=((-nx)/2:(nx)/2-1)/nx;            
y=((-ny)/2:(ny)/2-1)/nx; % Divide by nx (not ny) to ensure equal step size in x and y
[X,Y]=meshgrid(x,y);
Rad=sqrt(X.^2+Y.^2);
Angle=angle(X+1i.*Y)+pi; %Matrix with all Angles starting left-center

%% Normalize field
PreNormField=sum(sum(field.*conj(field)));
NormField=sqrt(real(PreNormField));
field=field/NormField;

%% Compare reconstructed field to signal without reference
figure(40)
subplot(121);
imagesc(abs(field).^2);
axis square
title("Intensity of reconstructed field")
subplot(122);
imagesc(signal_img);
title("Measured intensity")
axis square
colormap(gray);

%% Center, rotate, filter reconstructed field

% Center the field
field_centered = circshift(field,nx/2-cx_f,2);
field_centered = circshift(field_centered,ny/2-cy_f,1);
% figure(12)
% intensity_phase_plot(field_centered);

% Find center of mass of real space
if find_center_of_masses_flag
    [cx_i,cy_i] = center_of_mass(field_centered,center_of_mass_threshold);
    d_cx_i = cx_i-nx/2 % Difference to currently chosen center coordinate
    d_cy_i = cy_i-ny/2 % Difference to currently chosen center coordinate
end

% Rotate field
if rotate_field_flag
    field_centered = imrotate(field_centered,rotation_angle,'crop');
end

% Filter out a ring shape (reduce effect of stray light (not in RCF) on the
% OAM spectrum reconstruction)
field_centered(Rad<radi) = 0;
field_centered(Rad>rado) = 0;

% Visualize reconstructed, centered, rotated, filtered field
figure(13)
intensity_phase_plot(field_centered);
ax1 = gca;
title('Reconstructed field')

cbar = colorbar;
cbar.Label.String = "Phase";

% Show filtering ring
mask = (Rad> radi & Rad<radi+0.005);
mask(Rad>rado & Rad<rado+0.005) = 1;
ax2 = axes();
im_mask = imagesc(mask);
im_mask.AlphaData = im_mask.CData;
alpha scaled
ax2.ALim = [0 5];
xticks([])
yticks([])
ax2.Position = ax1.Position;
ax2.Color = 'None';
ax2.XColor = 'None';
ax2.YColor = 'None';
colormap(ax2, 'white')

%% OAM Spectrum

% Calculate OAM spectrum
OAM_lim = 30; % Consider OAM spectrum from -OAM_lim to OAM_lim
Spectrum = OAM_Spectrum(field_centered,Angle,Rad,OAM_lim); % Calculate OAM spectrum

% Clean up spectrum & Plot
phases = mod(-angle(Spectrum),2*pi);
phases = mod(phases-phases(OAM_lim+1), 2*pi); % zeroth component to 0 phase

% Correct linear phase due to rotation
if correct_linear_phase
    OAMs = -OAM_lim:OAM_lim;
    linear_phase = OAMs*rotation_angle/180*pi;
    phases = mod(phases+linear_phase,2*pi);
end

% Don't display phases for low intensity components
if filter_low_intensity_phases
    phases_filt = find(abs(Spectrum).^2<0.1*max(abs(Spectrum).^2));
    phases(phases_filt) = 0;
end

% Normalize spectrum to a max of 1
Spectrum = Spectrum./max(max(abs(Spectrum)));

% Plot spectrum
figure(14)
bar(-OAM_lim:OAM_lim,abs(Spectrum).^2,1)
if plot_phases
    hold on
    yyaxis right
    bar(-OAM_lim:OAM_lim,phases,0.5)
end
xlabel("OAM")
ylabel("Intensity (a.u.)")
ylim([0 1.3])
title('OAM Spectrum')
hold off
