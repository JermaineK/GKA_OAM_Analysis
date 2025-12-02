%% Constants & Flags

% Pixel size of images
nx = 320;
ny = 256;

% Max pixel value
norm = 4096;

% Whether to use special, wider sectors. Only implemented for datasets
% acquired using RCFs of length 7/12zT and 7/9 zT. Do not change here, it
% can be changed later when choosing the dataset
special_sectors_flag = 0;

%% Load data & specify settings

% SETTING EXPLANATIONS FOR DATASETS:
%
% Folder: Name of folder containing the measurement set
% bg: background image (no signal on camera)
% l: range of input OAM modes in measurement set
% n_sectors: Number of sorter sectors
% sector_angle_shift: Rotation of sorter sectors from default setting (radians)
% sector_ratio: The portion of the sector (0-1) considered in sorting.
%               Edges of the sectors are unsampled, if sector_ratio < 1.
%               This is only implemented for the special, wider sectors.
% x0, y0: Center coordinates of the field


% UNCOMMENT THE SELECTED DATASET SETTINGS FOR EVALUATION

% 2_3 zT
Folder = '2_3zT';
bg = double(imread(strcat(Folder,'\bg.png')))./norm;
l = -19:-12;
n_sectors = 3;
sector_angle_shift = -2*pi/3+pi/2-0.40;
x0 = 151.5639;
y0 = 126.1466;


% % 3_4 zT
% Folder = '3_4zT';
% bg = double(imread(strcat(Folder,'\background.png')))./norm;
% l = 14:21;
% n_sectors = 4;
% sector_angle_shift = .785;
% x0 = 151.5;
% y0 = 133.5;


% % 7_12 zT
% Folder = '7_12zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% l = -21:-10;
% n_sectors = 12;
% sector_angle_shift = -0.11 + pi/6;
% x0 = 155;
% y0 = 129;


% % 7_12 zT SPECIAL SECTORS
% Folder = '7_12zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% l = -18:-14;
% n_sectors = 5;
% special_sectors_flag = 1;
% sector_angle_shift = -1.45-5/12*2*pi + 0.06;
% sector_ratio = 1; % How big area of sector to consider (0-1)
% disp(strcat("Sector_ratio: ", num2str(sector_ratio)));
% x0 = 155;
% y0 = 129;


% % 7_9 zT
% Folder = '7_9zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% l = -21:-13;
% n_sectors = 9;
% sector_angle_shift = 2*pi/9*1.11;
% x0 = 165;
% y0 = 126.6189;


% % 7_9 zT SPECIAL SECTORS
% Folder = '7_9zT';
% bg = double(imread(strcat(Folder,'\bg.png')))./norm;
% l = -17:-14;
% special_sectors_flag = 1;
% n_sectors = 4;
% sector_angle_shift = 3*pi/9+0.84;
% sector_ratio = 1; % How big area of sector to consider (0-1)
% disp(strcat("Sector_ratio: ", num2str(sector_ratio)));
% x0 = 165;
% y0 = 126.61897;


% Load images
images = zeros(ny,nx,length(l));
for it = 1:length(l)
    images(:,:,it) = double(imread(strcat(Folder,'\l', num2str(l(it)), '.png')))./norm - bg;
end

% Shift the images of this data set a bit. (Use circshift to ensure the
% background noise is unaltered, reusing the same pixel values)
if strcmp(Folder,"7_9zT")
    images = circshift(images, [20 0 0]);
end

% Set negative values to 0 (If background reduction leads to negative values)?
images(images<0) = 0;

%% Linearize camera pixel response

% Pixel response linearization
c = [0.1562, 0.7375, -0.9792, 0.3774, 0.264, -0.1788].*1e-3;
images_original = images;
images = pix_linearize(c,images);

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

%% Generate crosstalk matrix

% Crosstalk matrix by summing the images multiplied by the binary sector masks
Crosstalk = zeros(n_sectors,length(l));
for it = 1:length(l)
    for jt = 1:n_sectors
        Crosstalk(jt,it) = sum(sum(images(:,:,it).*sector_masks(:,:,jt)));
    end
end

Crosstalk_raw = Crosstalk;

% Normalize crosstalk matrix (all cols sum up to 1)
% for it = 1:length(l)
%     Crosstalk(:,it) = Crosstalk(:,it)./sum(Crosstalk(:,it));
% end

% Normalize crosstalk matrix (all cols max is 1)
for it = 1:length(l)
    Crosstalk(:,it) = Crosstalk(:,it)./max(Crosstalk(:,it));
end

% Ratio of power localized in one sector (max sector/all sectors)
power_localization_ratio = zeros(1,length(l));
for it = 1:length(l)
    power_localization_ratio(it) = max(Crosstalk(:,it))./sum(Crosstalk(:,it));
end

Sorting_accuracy = mean(power_localization_ratio);
Sorting_accuracy_error = std(power_localization_ratio);

disp(strcat("Mean power localization ratio: ", num2str(Sorting_accuracy)))
disp(strcat("std of localization ratio: ", num2str(Sorting_accuracy_error)))

% Mean abs power localized to sectors (to track losses due to decreasing sector size)
disp(strcat("Mean abs power on diagonal: ", num2str(mean(max(Crosstalk_raw)))));
disp(strcat("Mean total abs power for a single mode: ", num2str(mean(sum(Crosstalk_raw)))));


%% Generate binary masks to plot sector lines
sector_outlines = zeros(ny,nx);
outline_width = 0.02;

for it = 1:n_sectors
    if ~special_sectors_flag % Standard sectoring
        ang = (it-1)*2*pi/n_sectors; % Angular position for line
        sector_outlines(Angle_shifted>=ang & Angle_shifted<ang + outline_width./(Rad*0.005)) = 1;
    else % Special sectoring (7/12zT)
        if strcmp(Folder,'7_12zT new\Xenics') % 7/12 sorting 5 sectors
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
        elseif strcmp(Folder, '7_9zT new\Xenics') % 7/12 sorting 4 sectors
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
if strcmp(Folder,"2_3zT") || strcmp(Folder,"3_4zT")
    n_plots_hor = 2;
    n_plots_ver = 4;
elseif strcmp(Folder,"7_12zT")
    n_plots_hor = 3;
    n_plots_ver = 4;
    if special_sectors_flag
        n_plots_hor = 1;
        n_plots_ver = 5;
    end
elseif strcmp(Folder,"7_9zT")
    n_plots_hor = 3;
    n_plots_ver = 3;
    if special_sectors_flag
        n_plots_hor = 1;
        n_plots_ver = 4;
    end
end

aspect_ratio = 1.6;
fig_width_cm = 8.8;
marg_left = 0.1;
marg_bot = 0.075;
offset_h = 0.01;
offset_v = 0.01;
scl = 0.1;
font_size = 9;

fig_width_cm = fig_width_cm*1.6;

fig = figure('units', 'centimeters', 'position', [30,10,fig_width_cm,fig_width_cm/aspect_ratio]);
set(fig, 'InvertHardcopy', 'off')
set(fig, 'Color', 'w')

% Colormap from white to color
col = [0 .482 0.671];
cmap_w_to_col = [linspace(1,col(1),256)',linspace(1,col(2),256)',linspace(1,col(3),256)'];

axes_vec = [];

it = 0;
% Plot fields with sector lines on top
for y = n_plots_ver-1:-1:0
    for x = 0:n_plots_hor-1
        it = it + 1;

        % Prepare axes
        % Positions of axes assuming square shaped, evenly spaced axes
        x_pos = marg_left + x*scl + x*offset_h;
        y_pos = (marg_bot + y*scl + y*offset_v)*aspect_ratio;
        sz_h = scl;
        sz_v = scl*aspect_ratio;
        ax = axes('Position', [x_pos, y_pos, sz_h, sz_v]);
        axes_vec = [axes_vec, ax];
        
        % Plot field
        img = images(:,:,it);
        square_img = img(:,floor(x0-ny/2)+1:floor(x0+ny/2));
        imagesc(square_img)
        colormap(ax,gray)
        xticks([])
        yticks([])
        ax.FontSize = font_size;

        % Plot sector lines
        ax2 = axes('Position', ax.Position);
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
        ax2.FontSize = font_size;

        xticks([])
        yticks([])
    end
end

% Colorbar for field intensity
cbar_img = colorbar(ax, 'Position', [marg_left-offset_h-sz_h/5, marg_bot*aspect_ratio, sz_h/5, ((n_plots_ver)*sz_v+(n_plots_ver-1)*offset_v*aspect_ratio)]);
cbar_img.Label.String = 'I (a.u.)';
cbar_img.Label.Position = [-.5, cbar_img.Limits(2)/2];
cbar_img.FontSize = font_size;
cbar_img.Ticks = cbar_img.Limits;
cbar_img.TickLabels = [0 1];


% Plot crosstalk

offset_cm = 0.05;

x_pos = marg_left + n_plots_hor*scl + n_plots_hor*offset_h + 0.1;
y_pos = (marg_bot + offset_cm)*aspect_ratio;
sz_v_cm = (n_plots_ver*sz_v+((n_plots_ver-1)*offset_v-offset_cm)*aspect_ratio);
sz_h_cm = sz_v_cm/size(Crosstalk,2)*size(Crosstalk,1)/aspect_ratio;
ax = axes('Position', [x_pos, y_pos, sz_h_cm, sz_v_cm]);

% Plot 
imagesc(Crosstalk')
yticks(1:length(l));
yticklabels(l);
ylabel('Input OAM')
xlabel('Sector')
ax.FontSize = font_size;
colormap(ax, cmap_w_to_col)

% Crosstalk colorbar
cbar_cm = colorbar('Position',[x_pos + sz_h_cm + offset_h, y_pos, sz_h/5, sz_v_cm]);
cbar_cm.Label.String = 'Power (a.u.)';
cbar_cm.Label.Position = [1, cbar_cm.Limits(2)/2];
cbar_cm.TickLabels = [0 1];
cbar_cm.FontSize = font_size;
cbar_cm.Ticks = cbar_cm.Limits;