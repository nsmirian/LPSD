%% take image from abr. beam camera in FL2 
% 
%  make sure that: (1) camera is not saturated
%                  (2) 
%  v1: swesch 26th Nov 2018: (1) TEST!!!! start to discriminate between SCR14FLFCOMP
%  v2: only comments: change to Christopher get_roi()
%  v3: add ta/ke image without streak to calculate time resolution: line 215
%  v4: ????
%  v5: 1st Nov

% 
% clear all
% close all
 
%    
function get_calibratedimage_OTR9FL2XTDS(block_laser, num_bgr, num_sig, msgLoc, comment)   
% comment         = 'FL2 PolariX measurement';
flag_save       = 1;
flag_act        = 1;

%
% num_bgr         = 10;  % number of single background images per set point
% num_sig         = 10;  % number of single measurement per set point

%
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');

%%

    % select cameras
    name_cam        = 'OTR9FL2XTDS'; 
    addr_cam        = ['FLASH.DIAG/CAMERA/', name_cam, '/'];
    name_script     = 'get_calibratedimage_OTR9FL2XTDS';
    fontSize        = 14;
    delay = doocsread([addr_cam, 'TRIGGERDELAYABS']);
    
    % polariX on off address
    addr_xtds_onoff = 'FLASH.DIAG/TIMINGINFO/FLFXTDS/ON_BEAM';

%% path

    addpath(genpath(pwd));

    path_tools          = '/home/ttflinac/user/swesch/Tools';
    addpath(genpath(path_tools));
    col                 = gen_ColorDefinitions;
    load('myColorMap.mat') % cmap


    % xfel tools
    path_xfel       = '/home/xfeloper/released_software/matlab/hlc_toolbox_common/';
    addpath(genpath(path_xfel));

    % statistic
    path_stats       = '/home/swesch/FFWD';
    addpath(path_stats);
    
    % get unstreaked image to calculate time resolution
    path_img         = '/home/ttflinac/user/swesch/FFWD/screen_image';
    addpath(path_img);
    
    % util_...
    path_lola       = '/home/ttflinac/doocs/source/matlab/LOLA-TDS/UltraShortPulse/';
    addpath(genpath(path_lola));
    load('myColorMap.mat') % cmap



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
       % poly parameter has to be switch on !!!!
    end

    % set image formats
    if (flag_act)
%         ddd_write = doocswrite([addr_cam, 'PIXELFORMAT.NUM'], 2);
%         ddd_write = doocswrite([addr_cam, 'FORMAT.OUT'], 1);
    end


%% prepare block laser for FLASH2

    % which laser
    name_laser          = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/LASER_SELECT.2'), 'data');
    addr_laser_block    = ['FLASH.DIAG/LASER.CONTROL/LASER', num2str(name_laser), '/BLOCK_LASER'];
    % addr_laser_block    = ['FLASH.DIAG/BEAMLINES/FLASH/LASER/BLOCK_LASER.FLASH3'];
    addr_laser_block    = ['FLASH.DIAG/BEAMLINES/FLASH/BLOCK_LASER.FLASH2'];
   
    % rep rate 
    rep_rate_macro      = getfield(doocsread('TTF2.UTIL/MAIN_PARAMETER/MACRO.REPRATE/VALUE'), 'data');
    bits_event7         = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/EVENT7'), 'data'); % FLASH2/FLASH3 '116'
    dividerA_event7     = bits_event7(4);
    rep_rate            = rep_rate_macro/(dividerA_event7+1);

%% prepare data structure

    % calibration data
    load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'time_res')
    load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')

    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    length_x            = size(ddd_read.data.val_val, 2);
    length_y            = size(ddd_read.data.val_val, 1);
    
    calib_x             = timecal_fspixel;
    time_axis           = timecal_fspixel * ((1:length_x)-length_x/2);
    
    calib_y             = ergcal; 
    erg_axis            = ergcal * ((1:length_y)-length_y/2);

    img_bgr             = zeros([size(ddd_read.data.val_val, 1), size(ddd_read.data.val_val, 2), num_bgr]);
    img_sig             = zeros([size(ddd_read.data.val_val, 1), size(ddd_read.data.val_val, 2), num_sig]);
    
    timeres_calib       = time_res;

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
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_bgr(:,:,jj)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_bgr(:,:,jj) == img_bgr(:,:,jj-1)
            ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            img_bgr(:,:,jj)       = ddd_read.data.val_val;
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
    img_bgr_mean = squeeze(mean(img_bgr, 3));
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
charge_7FL2XTDS     = zeros([1, num_sig]);
phase_TDS           = getfield(doocsread('FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE'),'data');
for jj = 1:num_sig

    %%% read spectrum x
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_sig(:,:,jj)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_sig(:,:,jj) == img_sig(:,:,jj-1)
            ddd_read        = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            img_sig(:,:,jj) = ddd_read.data.val_val;   
            display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(1/rep_rate)
        end
    end      
    
    % read charge in nC
    ddd_read = doocsread('FLASH.DIAG/TOROID/7FL2XTDS/CHARGE.FLASH2');
    charge_7FL2XTDS(jj) = ddd_read.data;

    pause(2/rep_rate)

end
    
display_message(msgLoc, ' - data taken');

% compute mean()
charge_7FL2XTDS_mean = mean(charge_7FL2XTDS); 

% clear stuff
clear tmp ddd_write ddd_read


%% save 

if (flag_save)
        save(['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp, '.mat'], ...
             'comment', 'timestamp', 'name_cam', 'img_bgr', 'img_sig', 'charge_7FL2XTDS', ...
             'timecal_fspixel', 'timecal_fspixel_err', 'timeres_calib', 'ergcal', 'ergcal_err', 'phase_TDS');
        display_message(msgLoc, [' - data saved']);
end


%% 

% analysis
img_filt    = img_sig;

num_sig     = size(img_sig, 3);
length_y    = size(img_sig, 1);
length_x    = size(img_sig, 2);
img_bgr_mean = squeeze(mean(img_bgr, 3));

x_com       = zeros([1, num_sig]);              y_com       = zeros([1, num_sig]);
x_var       = zeros([1, num_sig]);              y_var       = zeros([1, num_sig]);
x_fwhm      = zeros([1, num_sig]);              y_fwhm      = zeros([1, num_sig]);
x_profile   = zeros([num_sig, length_x]);       y_profile   = zeros([num_sig, length_y]);

for jj = 1:num_sig
    
%    tmp_img             = hlc_clean_image( squeeze(img_sig(:,:,jj)) - 0*img_bgr_mean );
    tmp_img             = get_ROI(squeeze(img_sig(:,:,jj)) - img_bgr_mean);
    img_filt(:,:,jj)    = tmp_img;
    
    % time
    tmp_profile         = mean(tmp_img);
%     tmp_profile(2010)=tmp_profile(2009); % remove hot pixel
%     tmp_profile         = medfilt1(tmp_profile, 1);
    [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);
    
    % erg
    tmp_profile         = mean(tmp_img, 2);
%     tmp_profile         = medfilt1(tmp_profile, 1);
    [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);
    
end


%% plot 1

fontSize        = 14;

erg_pos_good =  find(y_profile(end,:)>0); erg_pos_good = min(erg_pos_good):max(erg_pos_good);
time_pos_good = find(x_profile(end,:)>0); time_pos_good = min(time_pos_good):max(time_pos_good);



tmp2_img = squeeze(transpose(img_filt(:,:,end)));

% disp(' - ')
figure(1)
set(gcf, 'Position', [400, 1, 700, 700], 'Colormap', cmapZeroCubic)

    subplot(2,2,1)
        imagesc(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            erg_axis(erg_pos_good), flipud(tmp2_img(time_pos_good,erg_pos_good)') )
        grid on
        title(['single image ', name_cam], 'FontSize', fontSize, 'Interpreter', 'none')
        xlabel('time t /fs)', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}' , 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize,'YDir', 'normal')
        lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
        
    subplot(2,2,2)       
        p1 = plot(y_profile(end,erg_pos_good), -erg_axis(erg_pos_good) );
        grid on
        title(timestamp, 'FontSize', fontSize)
        xlabel('norm. erg density', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}', 'FontSize', fontSize)
        legend([p1(1)], {['rms: ' num2str(y_var(end), '%5.2f'), ' \cdot 10^{-3}', 10, ...
                'FWHM: ', num2str(y_fwhm(end), '%5.2f'), ' \cdot 10^{-3}']}, ...
                'FontSize', fontSize-2)        
        set(gca, 'FontSize', fontSize, 'YDir', 'normal');%, 'Ylim', -lim_y)        
        
    subplot(2,2,3)
    
    
 
    
    % just for quick fix jzemella
%     mu=sum( 1e-3*charge_7FLFMAFF_mean*1e-9*1e15 * x_profile(end,time_pos_good) .* time_axis(time_pos_good) )/sum(1e-3*charge_7FLFDUMP_mean*1e-9*1e15 * x_profile(end,time_pos_good));
%     var=sum(1e-3*charge_7FLFMAFF_mean*1e-9*1e15 * x_profile(end,time_pos_good) .* (time_axis(time_pos_good)-mu).^2 )/sum(1e-3*charge_7FLFDUMP_mean*1e-9*1e15 * x_profile(end,time_pos_good));
%     x_var=sqrt(var);
    
        p1 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            1e-3*charge_7FL2XTDS_mean*1e-9*1e15 * x_profile(end,time_pos_good));
        grid on
        title(comment, 'Interpreter', 'none')
        xlabel('time t /fs', 'FontSize', fontSize)        
        ylabel('peak current /kA', 'FontSize', fontSize)
%         legend([p1], {['rms mean: ' num2str(mean(x_var), '%5.2f'), ' fs', 10, ...
%                 'FWHM mean: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
%                 'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
%                 'FontSize', fontSize-2)      
        legend([p1], {['rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
                        'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
                        'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
                        'FontSize', fontSize-2)      
        set(gca, 'FontSize', fontSize, 'Xlim', lim_x)           
                %'FWHM: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
        %         

    subplot(2,2,4)
    MarkerSizz = 12;
    color1 = [0,0.4470, 0.7410];
    color2 = [0.85, 0.325, 0.098];
    yyaxis left
        p1 = plot(1:num_sig, x_var, 'kx', 'MarkerSize', MarkerSizz);
        hold on
        p2 = plot(1:num_sig, x_fwhm, 'x', 'Color', color1, 'MarkerSize', MarkerSizz);
        hold off
        ylabel('bunch length \sigma_t /fs', 'FontSize', fontSize)
        yyaxis right
           p3 = plot(1:num_sig, y_var, 'x', 'Color', color2, 'MarkerSize', MarkerSizz);
           ylabel('rms energy spread \sigma_E /10^{-3}', 'FontSize', fontSize)
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = color2;
        grid on
        xlabel('#measurement', 'FontSize', fontSize)
        
        legend([p1(1), p2(1), p3(1)], {['rms bu length, mean: ', num2str(mean(x_var), '%5.2f'), '\pm', num2str(std(x_var), '%5.2f') ' fs'], ...
                                       ['fwhm bu length, mean: ', num2str(mean(x_fwhm), '%5.2f'), '\pm', num2str(std(x_fwhm), '%5.2f') ' fs'], ...
                               ['rms E spread, mean: ', num2str(mean(y_var), '%5.2f'), '\pm', num2str(std(y_var), '%5.2f') ' \cdot 10^{-3}']}, ...
                          'FontSize', fontSize-2)
%                       legend([p1], {['rms mean: ' num2str(mean(x_var), '%5.2f'), ' fs', 10, ...
%                 'FWHM mean: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
%                 'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
%                 'FontSize', fontSize-2)      
%         legend([p1(1), p2(1)], {['x in ', calib_x_unit], ...
%                                ['y in ', calib_y_unit]}, ...
%                           'FontSize', fontSize-2)
        set(gca, 'FontSize', fontSize) 



    
     % add elog printing button
     filename=['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp, '.mat']
     text=[comment '\r\n', ...
             'data saved in:', filename]

     comment='data saved at '
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX LPSD@OTR9FL2XTDS', text));

end
