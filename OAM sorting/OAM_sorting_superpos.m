%% Constants & Parameters

nx = 320;
ny = 256;

norm = 4096;
%% Load data & specify settings

% SETTING EXPLANATIONS FOR DATASETS:
%
% Folder: Name of folder containing the measurement set
% bg: background image (no signal on camera)
% n_sectors: Number of sorter sectors
% sector_angle_shift: Rotation of sorter sectors from default setting (radians)
% sector_ratio: The portion of the sector (0-1) considered in sorting.
%               Edges of the sectors are unsampled, if sector_ratio < 1.
%               This is only implemented for the special, wider sectors.
% x0, y0: Center coordinates of the field

% UNCOMMENT THE SELECTED DATASET SETTINGS FOR EVALUATION


% 7_12 zT SPECIAL SECTORS. INPUT: L14, L15, L16, AMPLITUDES: 1, sqrt(2), sqrt(3)
Folder = '7_12zT';
bg = double(imread(strcat(Folder,'\bg.png')))./norm;
n_sectors = 5;
special_sectors_flag = 1;
sector_angle_shift = 6*pi/5-1.45-0.15; % 6 sectors, 5 modes
sector_ratio = 1; % How big area of sector to consider (0-1)
disp(strcat("Sector_ratio: ", num2str(sector_ratio)));
x0 = 155;
y0 = 129;
image = double(imread(strcat(Folder,'\Superpos\l-14_-15_-16_1_1.414_1.732.png')))./norm - bg;
spectrum_in = [0 3 1 0 2];
lvec = [-18 -16 -14 -17 -15];


% % 7_9 zT SPECIAL SECTORS. INPUT: L-15, L-16, L-17, AMPLITUDES: 1, 1, 1
% Folder = '7_9zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% special_sectors_flag = 1;
% n_sectors = 4;
% sector_angle_shift = 3*pi/9+0.84;
% sector_ratio = 1; % How big area of sector to consider (0-1)
% disp(strcat("Sector_ratio: ", num2str(sector_ratio)));
% x0 = 165;
% y0 = 126.61897;
% image = double(imread(strcat(Folder,'\Superpos\l-15_-16_-17_even.png')))./norm - bg;
% spectrum_in = [0 1 1 1];
% lvec = [-14 -15 -16 -17];


% % 7_9 zT SPECIAL SECTORS. INPUT: L-15, L-16, L-17, AMPLITUDES: 1, sqrt(2), sqrt(3)
% Folder = '7_9zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% special_sectors_flag = 1;
% n_sectors = 4;
% sector_angle_shift = 3*pi/9+0.84;
% sector_ratio = 1; % How big area of sector to consider (0-1)
% disp(strcat("Sector_ratio: ", num2str(sector_ratio)));
% x0 = 165;
% y0 = 126.61897;
% image = double(imread(strcat(Folder,'\Superpos\l-15_-16_-17_1_1.414_1.732.png')))./norm - bg;
% spectrum_in = [0 1 2 3];
% lvec = [-14 -15 -16 -17];



% Normalize spectrum_in
spectrum_in = spectrum_in./sum(spectrum_in);

% Set negative values to 0 (If background reduction leads to negative values)?
image(image<0) = 0;

% Shift the images of this data set a bit. (Use circshift to ensure the
% background noise is unaltered, reusing the same pixel values)
if strcmp(Folder,"7_9zT")
    image = circshift(image, [20 0 0]);
end

%% Linearize camera pixel response

% Pixel response linearization
c = [0.1562, 0.7375, -0.9792, 0.3774, 0.264, -0.1788].*1e-3;
image_original = image;
image = pix_linearize(c,image);

%% Coordinate system

% Generate coordinate system w.r.t center of the field
x = (1:nx)-x0;
y = (1:ny)-y0;
[X,Y] = meshgrid(x,y);
Rad=sqrt(X.^2+Y.^2);
Angle=angle(X+1i.*Y)+pi;

% Shifted Angle matrix (shift defined in load data section)
Angle_shifted = mod(Angle+sector_angle_shift,2*pi);
Angle_shifted = 2*pi-Angle_shifted; % Reverse the orientation (used in sector mask generation)


%% Binary Sector masks for sampling the field
sector_masks = zeros(ny,nx,n_sectors);
sector_angle = 2*pi/n_sectors;

% Angular masks
for it = 1:n_sectors
    % Standard
    mask = zeros(ny,nx);
    mask(Angle_shifted>(it-1)*sector_angle) = 1;
    mask(Angle_shifted>=it*sector_angle) = 0;
    sector_masks(:,:,it) = mask;
    
    % Replace the standard with special sectoring if special_sectors_flag = 1
    if special_sectors_flag % 7/12 special 5 sectors
        if strcmp(Folder,'7_12zT') % 7/12 sorting 5 sectors

            sector_spacing = 2*pi/6; % Instead of 2pi/5;
            sector_angle = sector_spacing*sector_ratio; % Smaller sectors
            dead_zone = (sector_spacing-sector_angle)/2;
            mask = zeros(ny,nx);
            if it < 4 % First 3 sectors as normal
                mask(Angle_shifted>(it-1)*sector_spacing+dead_zone) = 1;
                mask(Angle_shifted>=it*sector_spacing-dead_zone) = 0;
                sector_masks(:,:,it) = mask;
            else % Last 2 sectors shifted by additional 1/2 sector
                mask(Angle_shifted>((it+1/2)-1)*sector_spacing+dead_zone) = 1;
                mask(Angle_shifted>=(it+1/2)*sector_spacing-dead_zone) = 0;
                sector_masks(:,:,it) = mask;
            end

        elseif strcmp(Folder, '7_9zT') % 7/12 sorting 4 sectors

            sector_spacing = 4*pi/9;
            sector_angle = sector_spacing*sector_ratio; % Smaller sectors if sector_ratio < 1
            dead_zone = (sector_spacing-sector_angle)/2;
            mask = zeros(ny,nx);
            mask(Angle_shifted>(it-1)*sector_spacing+dead_zone) = 1;
            mask(Angle_shifted>=it*sector_spacing-dead_zone) = 0;
            sector_masks(:,:,it) = mask;
        end
    end
end


% Radial mask
Rmin = 82;
Rmax = 105;
Rmask = zeros(ny,nx);
Rmask(Rad>Rmin) = 1;
Rmask(Rad>Rmax) = 0;

sector_masks = sector_masks.*Rmask; % Multiply angular and radial masks

%% Calculate measured OAM spectrum

spectrum = zeros(1,n_sectors);
for it = 1:n_sectors
    spectrum(it) = sum(sum(image.*sector_masks(:,:,it)));
end

% Normalize spectrum
spectrum = spectrum./sum(spectrum);

% Error
mse = sum((spectrum_in-spectrum).^2)/n_sectors;
rmse = sqrt(mse);

%% Generate binary masks to plot sector lines
sector_outlines = zeros(ny,nx);
outline_width = 0.02;

for it = 1:n_sectors
    if ~special_sectors_flag % Standard sectoring
        ang = (it-1)*2*pi/n_sectors; % Angular position for line
        sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
    else % Special sectoring (7/12zT)
        if strcmp(Folder,'7_12zT') % 7/12 sorting 5 sectors
            if it < 4
                ang = (it-1)*2*pi/6; % Angular position for line
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
            elseif it == 4
                ang = (it-1)*2*pi/6; % End line for last standard sector
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
                ang = (it-1+1/2)*2*pi/6; % Start line for current sector
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
            elseif it == 5
                ang = (it-1+1/2)*2*pi/6; % Start line for current sector
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
                ang = (it+1/2)*2*pi/6; % End line for current sector
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
            end
        elseif strcmp(Folder, '7_9zT') % 7/12 sorting 4 sectors
            ang = (it-1)*4*pi/9;
            sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
            if it == 4 % End line
                ang = it*4*pi/9;
                sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
            end
        end
    end
end

sector_outlines(Rad<Rmin) = 0;
sector_outlines(Rad>=Rmin - 200*outline_width & Rad<Rmin) = 1;
sector_outlines(Rad>=Rmax & Rad<Rmax + 200*outline_width) = 1;
sector_outlines(Rad>Rmax + 200*outline_width) = 0;

%% Plotting

close all

% Define figure structure
aspect_ratio = 1.6;
font_size = 7;
fig_width_cm = 4.4;
fig_width_cm = fig_width_cm*1.6;
scl = 0.1875;

fig = figure('units', 'centimeters', 'position', [40,10,fig_width_cm,fig_width_cm/aspect_ratio]);
set(fig, 'InvertHardcopy', 'off')
set(fig, 'Color', 'w')

% Colormap from white to color
col = [0 .482 0.671];
cmap_w_to_col = [linspace(1,col(1),256)',linspace(1,col(2),256)',linspace(1,col(3),256)'];


% Camera image
ax1 = axes('Position', [0.1, 0.25, scl, scl*aspect_ratio]);
square_image = image(:,floor(x0-ny/2)+1:floor(x0+ny/2));
imagesc(square_image)
colormap(ax1,gray)
xticks([])
yticks([])
ax1.FontSize = font_size;


% Sector overlay
ax2 = axes('Position', ax1.Position);
im2 = imagesc(sector_outlines(:,floor(x0-ny/2)+1:floor(x0+ny/2)));
im2.AlphaData = im2.CData;
alpha scaled
ax2.ALim = [0 1];
colormap(ax2,col)
xticks([])
yticks([])
ax2.Color = 'None';
ax2.XColor = 'None';
ax2.YColor = 'None';


% Retrieved spectrum
ax3 = axes('Position', [0.5 0.25 scl*2 scl*aspect_ratio]);
[sorted_lvec, sorted_indices] = sort(lvec); % Sort OAM indices
col2 = [0.6706    0.1882         0];
b = bar(sorted_lvec, [spectrum_in(sorted_indices)', spectrum(sorted_indices)']);
b(1).FaceColor = col;
xticks(sorted_lvec)
box on
xlim([min(lvec)-0.6, max(lvec)+0.6])
xlabel('OAM')
ylabel({'Intensity';'(a.u.)'})
ylim([0 max([spectrum_in, spectrum])*1.1])
text(ax3.XLim(1),ax3.YLim(2)*1.2, strcat("RMSE: ", num2str(rmse,2)),'FontSize',font_size) % Error text
ax3.FontSize = font_size;
