function sigma_x = get_image_fl2_function_v3(num_bgr, num_sig, block_laser, msgLoc)
%% take image from abr. beam camera in FLF 
% 
%  make sure that: (1) camera is not saturated
%                  (2) 
%  v1: swesch 26th Nov 2018: (1) TEST!!!! start to discriminate between SCR14FLFCOMP
%  v2: only comments: change to Christopher get_roi()
%  v3: mflorian: change from sigma_x to sigma_y x -> x is energy, y is time


% clear all
% close all

% use Gauss Fit
method_utilGauss = 1;

% to suppress an matlab error - I do not know why.
ones(10)*ones(10);
    
% select cameras
% addr_cam        = 'FLASH.DIAG/CAMERA/14FLFCOMP/';
% addr_cam        = 'FLASH.DIAG/CAMERA/8FLFMAFF/';
% addr_cam        = 'FLASH.DIAG/CAMERA/5FLFDIAG/';
addr_cam        = 'FLASH.DIAG/CAMERA/OTR9FL2XTDS/';
% addr_cam        = 'FLASH.DIAG/FLASHFWDCAM12.CAM/SCREEN_DS_WALL/';
% addr_cam        = 'FLASH.DIAG/FLASHFWDCAM12.CAM/TARGET_SCREEN/';
% addr_cam        = 'FLASH.DIAG/FLASHFWDCAM7.CAM/EspecHighres/';
% addr_cam        = 'FLASH.DIAG/FLASHFWDCAM7.CAM/3.5_FLFDIAG/';
%addr_cam        = 'FLASH.DIAG/FLASHFWDCAM12.CAM/TARGET_WALL/';
delay = doocsread([addr_cam, 'TRIGGERDELAYABS']);

comment         = 'OTR9FL2XTDS image to get time resolution.';
% 
flag_save       = 1;
flag_act        = 1;

%
% num_bgr         = 5;  % number of single background images per set point
% num_sig         = 10;  % number of single measurement per set point

% 	block_laser = 0; % 1 = use injector laser; 0 = use trigger delay

%% 
name_script     = 'Intrinsic beam size measurement FL2';
name_cam        = getfield(doocsread([addr_cam ,'NAME']), 'data');
fontSize        = 14;

% time
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');

% util_...
path_lola       = '/home/ttflinac/doocs/source/matlab/LOLA-TDS/UltraShortPulse/';
addpath(genpath(path_lola));
load('myColorMap.mat') % cmap


% xfel tools
path_xfel       = '/home/xfeloper/released_software/matlab/hlc_toolbox_common/';
addpath(genpath(path_xfel));

% statistic
path_stats       = '/home/swesch/FFWD';
addpath(path_stats);

display_message(msgLoc, [name_script, '(): start']);

%% prepare camera 

% switch on camera

% switch off ROIs
if (flag_act)
    ddd_write = doocswrite([addr_cam, 'ROI_SPECTRUM.ON'], 0);
    ddd_write = doocswrite([addr_cam, 'ROI2_SPECTRUM.ON'], 0);
end

% switch on spectrum
if (flag_act)
    ddd_write = doocswrite([addr_cam, 'SPECTRUM.ON'], 1);
end

ddd_read        = doocsread([addr_cam, 'X.POLY_SCALE']);
calib_x         = 1e-3 * ddd_read.data(3); % m/pixel
calib_x_plot    = 1e3; % mm
calib_x_unit    = 'mm';
if calib_x == 0
    calib_x         = 1;
    calib_x_plot    = 1;    
    calib_x_unit    = 'pixel';
end

ddd_read        = doocsread([addr_cam, 'Y.POLY_SCALE']);
calib_y         = 1e-3 * ddd_read.data(3); % m/pixel
calib_y_plot    = 1e3; % mm
calib_y_unit    = 'mm';
if calib_y == 0
    calib_y         = 1;
    calib_y_plot    = 1;
    calib_y_unit    = 'pixel';
end



%% prepare block laser for FLASH2

% which laser
name_laser          = getfield(doocsread(['FLASH.DIAG/TIMER/FLASHCPUTIME1.0/LASER_SELECT.2']), 'data');
% if (name_laser == 0)
%    display_message(msgLoc, ' - no LASER for FLASH3 selected ... please stop here ... \n')
% end
% addr_laser_block    = ['FLASH.DIAG/LASER.CONTROL/LASER', num2str(name_laser), '/BLOCK_LASER'];

 addr_laser_block    = ['FLASH.DIAG/BEAMLINES/FLASH/BLOCK_LASER.FLASH2'];

% macro rep rate 
rep_rate_macro      = getfield(doocsread('TTF2.UTIL/MAIN_PARAMETER/MACRO.REPRATE/VALUE'), 'data');
bits_event7         = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/EVENT7'), 'data'); % FLASH2/FLASH3 '116'
dividerA_event7     = bits_event7(4);
rep_rate            = rep_rate_macro/(dividerA_event7*9+1);


%% prepare data structure


ddd_read            = doocsread([addr_cam, 'IMAGE_EXT']);
length_x            = size(ddd_read.data.val_val, 2);
length_y            = size(ddd_read.data.val_val, 1);

img_bgr             = zeros([num_bgr, size(ddd_read.data.val_val, 1), size(ddd_read.data.val_val, 2)]);
img_sig             = zeros([num_sig, size(ddd_read.data.val_val, 1), size(ddd_read.data.val_val, 2)]);

display_message(msgLoc, ' - setup done');

%% take background 
    
% block laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 1);
        display_message(msgLoc, [' - block LASER', name_laser]);  
        pause(2)         
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 10000);
    display_message(msgLoc, [' - changed camera delay']);  
    pause(2)
end

% take bgr
for jj = 1:num_bgr 
    
    % read
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT']);
    img_bgr(jj,:,:)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_bgr(jj,:,:) == img_bgr(jj-1,:,:)
            ddd_read            = doocsread([addr_cam, 'IMAGE_EXT']);
            img_bgr(jj,:,:)       = ddd_read.data.val_val;
            display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(10/rep_rate)
        end
    end
    
    pause(1/rep_rate)
    
end
    
% compute mean()
if num_bgr == 1
    img_bgr_mean = squeeze(img_bgr);
else
    img_bgr_mean = squeeze(mean(img_bgr));
end


display_message(msgLoc, ' - background data taken');

%% take data


% unblock laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 0); % unblock laser 
        display_message(msgLoc, [name_script, '(): Unblock LASER', name_laser]);
        pause(0.5)       
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], delay.data);
    display_message(msgLoc, [' - changed camera delay back']);  
    pause(2)
end


    
  
% loop
for jj = 1:num_sig

    %%% read spectrum x
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT']);
    img_sig(jj,:,:)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_sig(jj,:,:) == img_sig(jj-1,:,:)   % it can be bettter, we should change to time stamp 
            ddd_read        = doocsread([addr_cam, 'IMAGE_EXT']);
            img_sig(jj,:,:) = ddd_read.data.val_val;   
            display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(1/rep_rate)
        end
    end        

    pause(2/rep_rate)

end
    
display_message(msgLoc, ' - data taken');


ddd_read = doocsread('FLASH.DIAG/TOROID/18FL2EXTR/CHARGE.FLASH2'); %I change to 2  it was 3
charge_18FLF2EXTR = ddd_read.data;

ddd_read = doocsread('FLASH.DIAG/TOROID/5FLFEXTR/CHARGE.FLASH2');
charge_5FLFEXTR = ddd_read.data;

ddd_read = doocsread('FLASH.DIAG/TOROID/7FLFMAFF/CHARGE.FLASH2');
charge_7FLFMAFF = ddd_read.data;



% clear stuff
clear tmp ddd_write ddd_read


%% save 

if (flag_save)
        save(['image_reso/image_', name_cam, '_', timestamp, '.mat']);
        display_message(msgLoc, [' - data saved']);
end

%% 

% analysis
img_filt    = img_sig;

x_com       = zeros([1, num_sig]);              y_com       = zeros([1, num_sig]);
x_var       = zeros([1, num_sig]);              y_var       = zeros([1, num_sig]);
x_fwhm      = zeros([1, num_sig]);              y_fwhm      = zeros([1, num_sig]);
x_profile   = zeros([num_sig, length_x]);       y_profile   = zeros([num_sig, length_y]);

for jj = 1:num_sig
    
    tmp_img             = hlc_clean_image( squeeze(img_sig(jj,:,:))-0.*img_bgr_mean );
    img_filt(jj,:,:)    = tmp_img;
    
    tmp_profile         = mean(tmp_img);
%     tmp_profile         = medfilt1(tmp_profile, 1);
%     [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);
    
    if method_utilGauss
        n       = size(tmp_profile,2); % number of bins
        a       = calib_x * ( (1:n) - n/2 ); % axis
        [par_gauss, ~] = util_gaussFit(a,tmp_profile, 1);
        x_com(jj) = par_gauss(2);
        x_var(jj) = par_gauss(3);
        [~, ~, x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);
    %     x_gauss = par_avg_gaussfit(1)*exp(-((x_axis-par_avg_gaussfit(2))/sqrt(2)/par_avg_gaussfit(3)).^2);
    else
        [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);
    end
    
    tmp_profile         = mean(tmp_img, 2);
%     tmp_profile         = medfilt1(tmp_profile, 1);
%     [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);
    if method_utilGauss
        n       = size(tmp_profile,1); % number of bins
        a       = calib_y * ( (1:n) - n/2 ); % axis
        [par_gauss, x_gauss] = util_gaussFit(a,tmp_profile', 1);
        y_com(jj) = par_gauss(2);
        y_var(jj) = par_gauss(3);
        [~, ~, y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);
    %     x_gauss = par_avg_gaussfit(1)*exp(-((x_axis-par_avg_gaussfit(2))/sqrt(2)/par_avg_gaussfit(3)).^2);
    else
        [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);
    end
    
end


%% plot 1

% disp(' - ')
figure(1)
set(gcf, 'Position', [800, 1, 800, 800], 'Colormap', cmapZeroCubic)

    subplot(2,2,1)
        imagesc(calib_x_plot*x_axis, calib_y_plot*y_axis, squeeze(img_filt(1,:,:)))
        grid on
        title(['single image ', name_cam], 'FontSize', fontSize, 'Interpreter', 'none')
        xlabel(['horizontal position x (', calib_x_unit, ')'], 'FontSize', fontSize)        
        ylabel(['vertical position y (', calib_y_unit, ')'], 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize)
        lim_y = get(gca, 'YLim');lim_x = get(gca, 'XLim');
        
    subplot(2,2,2)       
        p1 = plot(calib_y_plot*y_axis, 1/calib_y_plot*y_profile);
        grid on
        title(timestamp, 'FontSize', fontSize)
        ylabel(['norm. density \rho_y (', calib_y_unit, '^{-1})'], 'FontSize', fontSize)        
        xlabel(['vertical position y (', calib_y_unit, ')'], 'FontSize', fontSize)
        legend([p1(1)], {['rms: ' num2str(calib_y_plot*mean(y_var), '%5.2f'), ' ', calib_y_unit, 10, ...
                'FWHM: ', num2str(calib_y_plot*mean(y_fwhm), '%5.2f'), ' ', calib_y_unit]}, ...
                'FontSize', fontSize-2)        
        set(gca, 'FontSize', fontSize, 'YDir', 'normal', 'Xlim', lim_y)        
        
    subplot(2,2,3)
        p1 = plot(calib_x_plot*x_axis, 1/calib_x_plot*x_profile);
        grid on
        title(comment, 'Interpreter', 'none')
        xlabel(['horizontal position x (', calib_x_unit, ')'], 'FontSize', fontSize)        
        ylabel(['norm. density \rho_x (', calib_x_unit, '^{-1})'], 'FontSize', fontSize)
        legend([p1], {['rms: ' num2str(calib_x_plot*mean(x_var), '%5.2f'), ' ', calib_x_unit, 10, ...
                'FWHM: ', num2str(calib_x_plot*mean(x_fwhm), '%5.2f'), ' ', calib_x_unit]}, ...
                'FontSize', fontSize-2)        
        set(gca, 'FontSize', fontSize, 'Xlim', lim_x)           
        

    subplot(2,2,4)
        p1 = plot(1:num_sig, calib_x_plot*x_var, 'b.');
        hold on
           p2 = plot(1:num_sig, calib_y_plot*y_var, '.r');
        hold off
        grid on
        xlabel('measurement', 'FontSize', fontSize)
        ylabel('rms beamsize \sigma (mm)', 'FontSize', fontSize)
        legend([p1(1), p2(1)], {['x in ', calib_x_unit], ...
                               ['y in ', calib_y_unit]}, ...
                          'FontSize', fontSize-2)
        set(gca, 'FontSize', fontSize)
        
    sigma_x = calib_x_plot*mean(x_var);
    
    uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX image for time resolution', comment));

% %% plot 2
% 
% if strcmp(addr_cam, 'FLASH.DIAG/CAMERA/14FLFCOMP/')
%     design_R16 = 0.206; % in m
%     calib_x_unit = '/ 10^{3}';
%     calib_x_plot = 1000/design_R16; % in permille
%     label_x_axis = ['rel. energy \delta ', calib_x_unit];
%     label_y_axis = ['norm. density \rho_\delta ', calib_x_unit];
% else
%     design_R16 = 1; % in mm
%     label_x_axis = ['horizontal position x (', calib_x_unit, ')'];
%     label_y_axis = ['norm. density \rho_x (', calib_x_unit, '^{-1})'];
% end
% 
% axis_x      = calib_x_plot*x_axis;
% profile_x   = 1/calib_x_plot*x_profile(1,:);
% 
% 
% axis_y      = calib_y_plot*y_axis;
% profile_y   = 1/calib_y_plot*y_profile(1,:);
% 
% 
% 
% disp(' - ')
% figure(2)
% set(gcf, 'OuterPosition', [800, 1, 800, 800], 'Colormap', cmapZeroCubic)
% 
%     subplot(2,2,1)
%         imagesc(axis_x, axis_y, squeeze(img_filt(1,:,:)))
%         grid on
%         title(['single image ', name_cam], 'FontSize', fontSize, 'Interpreter', 'none')
%         xlabel(label_x_axis, 'FontSize', fontSize)        
%         ylabel(['vertical position y (', calib_y_unit, ')'], 'FontSize', fontSize)
%         set(gca, 'FontSize', fontSize)
%         lim_y = get(gca, 'YLim');lim_x = get(gca, 'XLim');
%         
%     subplot(2,2,2)       
%         p1 = plot(profile_y, axis_y);
%         grid on
%         title(timestamp, 'FontSize', fontSize)
%         xlabel(['norm. density \rho_y (', calib_y_unit, '^{-1})'], 'FontSize', fontSize)        
% %         ylabel(['vertical position y (', calib_y_unit, ')'], 'FontSize', fontSize)
%         legend([p1], {['rms: ' num2str(calib_y_plot*mean(y_var), '%5.2f'), ' ', calib_y_unit, 10, ...
%                 'FWHM: ', num2str(calib_y_plot*mean(y_fwhm), '%5.2f'), ' ', calib_y_unit]}, ...
%                 'FontSize', fontSize-2)  
%         set(gca, 'FontSize', fontSize, ...
%                  'YDir', 'reverse', ...
%                  'Ylim', lim_y, ...
%                  'Xlim', [0, 1.05*max(profile_y)])        
%         
%     subplot(2,2,3)
%         p1 = plot(axis_x, profile_x);
%         grid on
%         xlabel(label_x_axis, 'FontSize', fontSize)        
%         ylabel(label_y_axis, 'FontSize', fontSize)
%         legend([p1], {['rms: ' num2str(calib_x_plot*mean(x_var), '%5.2f'), ' ', calib_x_unit, 10, ...
%                 'FWHM: ', num2str(calib_x_plot*mean(x_fwhm), '%5.2f'), ' ', calib_x_unit]}, ...
%                 'FontSize', fontSize-2)        
%         set(gca, 'FontSize', fontSize, 'Xlim', lim_x)           
%    
%         
%         
%         
