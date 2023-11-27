%% script for PolariX time calibration for screen OTR9FL2XTDS
%
% swesch - 8th Sep 2019 - v1: first script 
% 
% v2: add timecal_streak to calculate time resolution, last 1st Nov

function get_timecalibration_OTR9FL2XTDS(start_actuator, end_actuator, num_sig)
    
%% config

%     comment         = 'FL2 PolariX';

    
    flag_pause      = 1; % make a pause between reach goal current and take measurement

% 	bgr_method = 0; % 1 = use injector laser; 0 = use trigger delay

    fontSize        = 14;
%% 

    % addresses
    name_script             = 'get_timecalibration_OTR9FL2XTDS';
    name_cam                = 'OTR9FL2XTDS';
    addr_cam                = ['FLASH.DIAG/CAMERA/', name_cam, '/'];
    delay = doocsread([addr_cam, 'TRIGGERDELAYABS']);
    
    pixelformat=doocswrite([addr_cam, 'FORMAT.OUT'], 0); %until IAS is used: set pixelformat out to mono8

    % address magnet directly
    name_actuator           = 'XTDS phase';
    addr_actuator_set       = 'FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
    addr_actuator_rbv       = 'FLASH.RF/LLRF.CONTROLLER/FORW.SLED.POLARIX/PHASE.SAMPLE';

    % time
    timestamp               = datestr(clock, 'yyyy-mm-ddTHHMMSS');
    
    % PolariX
    frequency_XTDS          = 11.9888e9; % in Hz  
    ddd_read                = doocsread('FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.AMPL');
    amplitude_XTDS          = ddd_read.data; % in %
    addr_xtds_onoff         = 'FLASH.DIAG/TIMINGINFO/FLFXTDS/ON_BEAM';


    % get unstreaked image to calculate time resolution
    path_img         = '/home/ttflinac/user/swesch/FFWD/screen_image';
    addpath(path_img);


%% path

    path_tools          = '/home/ttflinac/user/swesch/Tools';
    addpath(genpath(path_tools));
    col                 = gen_ColorDefinitions;
    load('myColorMap.mat') % cmap

%% prepare camera 

    % switch off ROIs
    if (flag_act)
        ddd_write = doocswrite([addr_cam, 'ROI_SPECTRUM.ON'], 0);
        ddd_write = doocswrite([addr_cam, 'ROI2_SPECTRUM.ON'], 0);
    end

    % switch off: BG subtraction
    if (flag_act)
        ddd_write = doocswrite([addr_cam, 'SUBSTR.ON'], 0);
    end    
    
    % switch on spectrum
    if (flag_act)
        ddd_write = doocswrite([addr_cam, 'SPECTRUM.ON'], 1);       
    end

    
    
%% prepare data structure

    ddd_read        = doocsread([addr_cam, 'IMAGE_EXT']);
    cam_spec        = ddd_read.data;

    % background (shot, spectrum)
    bgr_x           = zeros(num_bgr, size(cam_spec.val_val,2));
    bgr_y           = zeros(num_bgr, size(cam_spec.val_val,1));

    % signal (scan_point, shot, spectrum)
    raw_spec_x      = zeros(num_actuator, num_sig, size(cam_spec.val_val,2));
    corr_spec_x     = raw_spec_x;
    raw_spec_y      = zeros(num_actuator, num_sig, size(cam_spec.val_val,1));
    corr_spec_y     = raw_spec_y;
    

    % center of mass
    pos_CoM_x       = zeros(num_actuator, num_sig);
    pos_CoM_y       = zeros(num_actuator, num_sig);

    % scales 
    scale_x         = abs(cam_spec.scale_x); 
    scale_y         = abs(cam_spec.scale_y);
    if scale_x == 1
        scale_x         = 0.006184; % mm/pixel
    end
    if scale_y == 1
        scale_y         = 0.006352; % mm/pixel 
    end
%     scale_y         = 0.0109; % mm/pixel 
%     scale_x         = 0.0114; % mm/pixel


    % current rbv
    actuator_rbv    = zeros([1, num_actuator]);

    % charge
    charge_7FLFMAFF = zeros(num_actuator, num_sig);
    charge_7FLFDUMP = zeros(num_actuator, num_sig);

%% prepare block laser for FLASH2

    % which laser
    name_laser          = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/LASER_SELECT.2'), 'data');
    %addr_laser_block    = ['FLASH.DIAG/LASER.CONTROL/LASER', num2str(name_laser), '/BLOCK_LASER'];
    addr_laser_block    = ['FLASH.DIAG/BEAMLINES/FLASH/BLOCK_LASER.FLASH2'];
   
    % rep rate 
    rep_rate_macro      = getfield(doocsread('TTF2.UTIL/MAIN_PARAMETER/MACRO.REPRATE/VALUE'), 'data');
    bits_event7         = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/EVENT7'), 'data'); % FLASH2/FLASH3 '116'
    dividerA_event7     = bits_event7(4);
    rep_rate            = rep_rate_macro/(dividerA_event7+1);

%% prepare actuator scan list


    % reference setpoint
    ref_actuator_set    = 1e-3*round(1e3*getfield(doocsread(addr_actuator_set), 'data'));
    ref_actuator_rbv    = 1e-3*round(1e3*getfield(doocsread(addr_actuator_rbv), 'data'));

    % scan list
    scan_list_set       = linspace(start_actuator, end_actuator, num_actuator);
    % rbv container
    scan_list_rbv       = zeros(num_actuator, num_sig);


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
    
    % read x
    ddd_read            = doocsread([addr_cam, 'SPECTRUM.X.TD']);
    bgr_x(jj,:)         = ddd_read.data.d_gspect_array_val;
    
    % compared it with last measurement
    if jj > 1
        while bgr_x(jj,:) == bgr_x(jj-1,:)
            ddd_read        = doocsread([addr_cam, 'SPECTRUM.X.TD']);
            bgr_x(jj,:)     = ddd_read.data.d_gspect_array_val;
            display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(10/rep_rate)
        end
    end
    
    % read y 
    ddd_read            = doocsread([addr_cam, 'SPECTRUM.Y.TD']);
    bgr_y(jj,:)         = ddd_read.data.d_gspect_array_val;    
    
    % compared it with last measurement
    if jj > 1
        while bgr_y(jj,:) == bgr_y(jj-1,:)
            ddd_read        = doocsread([addr_cam, 'SPECTRUM.Y.TD']);
            bgr_y(jj,:)     = ddd_read.data.d_gspect_array_val;
            display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(10/rep_rate)
        end
    end    
    
    pause(1/rep_rate)
    
end
    
% compute mean()
bgr_spec_x_mean = mean(bgr_x);
bgr_spec_y_mean = mean(bgr_y);

% take one img
ddd_read        = doocsread([addr_cam, 'IMAGE_EXT']);
bgr_img         = ddd_read.data.val_val;

% plot
figure(1)
set(gcf, 'Position', [800, 1, 800, 700], 'Colormap', cmapZeroCubic)

    subplot(2,2,1)
        imagesc(bgr_img)
        title(['background image', name_cam], 'FontSize', fontSize)
        ylabel('vertical axis (pixel)', 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize)        

    subplot(2,2,2)
        plot(bgr_spec_y_mean, 1:length(bgr_spec_y_mean))
        grid on
        xlabel('projection y', 'FontSize', fontSize)
        title(timestamp, 'FontSize', fontSize)
        ylim([1, length(bgr_spec_y_mean)])
        set(gca, 'YDir', 'Reverse', 'FontSize', fontSize)
        
    subplot(2,2,3)
        plot(bgr_spec_x_mean)
        grid on
        xlabel('horizontal axis (pixel)', 'FontSize', fontSize)
        ylabel('projection x', 'FontSize', fontSize)
        xlim([1, length(bgr_spec_x_mean)])
        set(gca, 'FontSize', fontSize) 


display_message(msgLoc, ' - background data taken');


%% scan loop
figure(2)
display_message(msgLoc, [name_script, '(): Start scan ...']);

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

% turn XTDS on:
tmp      = doocswrite(addr_xtds_onoff, 1);

for ii = 1:num_actuator % scan points
    
       
    %%%%%%%%%% set dipole current %%%%%%%%%%
    if (flag_act)
                        
        % set value 
        ddd_write = doocswrite(addr_actuator_set, scan_list_set(ii));
        display_message(msgLoc, [name_script, '(): Set actuator ', name_actuator, ' to ' , num2str(scan_list_set(ii), '%5.3f')]);
        
        % wait for set = rbv
        display_message(msgLoc, [name_script, '(): Actuator ', name_actuator, ' set.']);
        
        if (flag_pause)
            pause(5)
        end
        
    end   
    
    
    %%%%%%%%%%  take data %%%%%%%%%%

    for jj = 1:num_sig
        
        %%% read actuator readback
        ddd_read                = doocsread(addr_actuator_rbv);
        scan_list_rbv(ii,jj)    = ddd_read.data;
        
        %%% read spectrum x
        ddd_read                = doocsread([addr_cam, 'SPECTRUM.X.TD']);
        raw_spec_x(ii,jj,:)     = ddd_read.data.d_gspect_array_val;
        % compared it with last measurement
        if jj > 1
            while raw_spec_x(ii,jj,:) == raw_spec_x(ii,jj-1,:)
                ddd_read            = doocsread([addr_cam, 'SPECTRUM.X.TD']);
                raw_spec_x(ii,jj,:) = ddd_read.data.d_gspect_array_val;   
                display_message(msgLoc, [name_script, '(): same data ... wait ...']); 
                pause(1/rep_rate)
            end
        end        
        
        % remove bg x
        tmp                     = squeeze(raw_spec_x(ii,jj,:));
        corr_spec_x(ii,jj,:)    = tmp' - bgr_spec_x_mean;
 
        % filter (mean), fit (asymmetric gauss) and position (from fit)
        tmp                     = medfilt1(corr_spec_x(ii,jj,:), 3);
        tmp                     = util_gaussFit(1:length(tmp), tmp, 1, 1);
        pos_CoM_x(ii,jj)        = tmp(2);
        
        
        
        %%% read spectrum Y
        ddd_read                = doocsread([addr_cam, 'SPECTRUM.Y.TD']);
        raw_spec_y(ii,jj,:)     = ddd_read.data.d_gspect_array_val;      
        
        % remove bg y
        tmp                     = squeeze(raw_spec_y(ii,jj,:));
        corr_spec_y(ii,jj,:)    = tmp' - bgr_spec_y_mean;
        
        % filter (mean), fit (asymmetric gauss) and position (from fit)
        tmp                     = medfilt1(corr_spec_y(ii,jj,:), 3);
        [tmp, tmp2]             = util_gaussFit(1:length(tmp), tmp, 1, 1);
        pos_CoM_y(ii,jj)        = tmp(2);        
        
        
        pause(2/rep_rate)
        
    end
    
    % plot
    subplot(2,1,1)
        plot(scan_list_set, pos_CoM_x, 'b.')
        title([name_script, '() - ', timestamp], 'Interpreter', 'none', 'FontSize', fontSize)
%         xlabel(addr_actuator_set, 'Interpreter', 'none')
        ylabel('beam position <x> (pixel)', 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize) 
                
    subplot(2,1,2)
        plot(scan_list_set, pos_CoM_y, 'b.')
%         title([name_script, '() - ', timestamp], 'Interpreter', 'none')
        xlabel(addr_actuator_set, 'Interpreter', 'none', 'FontSize', fontSize)
        ylabel('beam position <y> (pixel)', 'FontSize', fontSize)        
        set(gca, 'FontSize', fontSize) 

end
display_message(msgLoc, [name_script, '(): Scan ended.']);


% restore actuator reference value 
if (flag_act)
    tmp = mean([start_actuator, end_actuator]);
    ddd_write = doocswrite(addr_actuator_set, tmp);
    display_message(msgLoc, [name_script, '(): Set actuator ', name_actuator, ' to ' , num2str(tmp, '%5.3f')]);
end


%% clear stuff
clear tmp tmp2 ddd_write ddd_read 



%% fit 


    % phase change
    new_scan_list   = mean(scan_list_rbv, 2)-mean(scan_list_set);
    new_scan_list   = scan_list_set-mean(scan_list_set);

    % x fit
    pos_CoM_x_mean  =  mean(pos_CoM_x, 2);  % pixel
    pos_CoM_x_std   =  std(pos_CoM_x, 1, 2); 
    [par_x, yFit_x, parstd_x] = util_polyFit(new_scan_list, pos_CoM_x_mean, 1, pos_CoM_x_std);
    
    timecal_fspixel = 1e15 * deg2rad(1/par_x(1))/(2*pi*frequency_XTDS);
    timecal_fspixel_err = 1e15 * deg2rad(parstd_x(1)/par_x(1)^2)/(2*pi*frequency_XTDS);
    
    timecal_streak  = 1/(3e8/(1e-3*scale_x) * deg2rad(1/par_x(1))/(2*pi*frequency_XTDS));
    timecal_streak_err  = 3e8/(1e-3*scale_x) * deg2rad(parstd_x(1)/par_x(1)^2)/(2*pi*frequency_XTDS);
    
    % y fit
    pos_CoM_y_mean  =  mean(pos_CoM_y, 2);  % pixel
    pos_CoM_y_std   =  std(pos_CoM_y, 1, 2); 
    [par_y, yFit_y, parstd_y] = util_polyFit(new_scan_list, pos_CoM_y_mean, 1, pos_CoM_y_std);

    
%% now take an image without streaking to calculate the time resolution

% turn off tds:
%     tmp      = doocswrite([addr_cam, 'GAINRAW'], 0);
%     tmp      = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 504.4);
%     tmp      = doocswrite([addr_cam, 'EXPOSURETIMEABS'], 18);
    tmp      = doocswrite(addr_xtds_onoff, 0);
    sigma_x  = get_image_fl2_function_v3(num_bgr, num_sig, block_laser, msgLoc);
    time_res = abs((sigma_x * 1e-3 / 3e8 ) / timecal_streak *1e15) ; % in m
    tmp      = doocswrite(addr_xtds_onoff, 1);

    
%% final plot

figure(4)
set(gcf, 'OuterPosition', [100, 100, 1000, 500])

    subplot(1,2,1)

    p1 = plot(new_scan_list, pos_CoM_x, 'b.');
    hold on
        p2 = plot(new_scan_list, yFit_x, '-r');
    hold off
    grid on
    title(['time calibration ', name_cam], 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel([name_actuator, ' (deg)'], 'FontSize', fontSize)
    ylabel([name_cam, ' beam position <y> (pixel)'], 'FontSize', fontSize)
    legend([p1(1), p2], ...
            {'data', ['fit: ', num2str(1/par_x(1)*1e3, '%5.1f'), ' +/- ', num2str(parstd_x(1)/par_x(1)^2*1e3, '%5.1f'), ' mdeg/pixel']}, ...
            'Location', 'NorthWest', 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize) 
    
    subplot(1,2,2)
    p1 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*pos_CoM_x, 'b.');
    hold on
        p2 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*yFit_x, '-r');
    hold off
    grid on
    title(timestamp, 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel([name_actuator, ' (ps)'], 'FontSize', fontSize)
    ylabel('beam position <y> (mm)', 'FontSize', fontSize) 
    legend([p1(1), p2], ...
        {'data', ['fit: ', num2str(timecal_fspixel, '%5.2f'), ' +/- ', num2str(timecal_fspixel_err, '%5.2f'), ' fs/pixel', 10 ,...
                            num2str(timecal_fspixel/scale_x, '%5.2f'), ' +/- ', num2str(timecal_fspixel_err/scale_x, '%5.2f'), ' fs/mm', 10 ,...
                           '(streak ', num2str(timecal_streak, '%5.1f'), ')', 10, ...
                           'Time res: ' num2str((time_res), '%5.1f'), ' fs']}, ...
        'Location', 'NorthWest', 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize) 
    
     % add elog printing button
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX time calibration', comment));


     % reset camera
%     tmp      = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 0);
%     tmp      = doocswrite([addr_cam, 'EXPOSURETIMEABS'], 1000);


%% save 

if (flag_save)
        save(['/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/', name_script(5:end), '_', timestamp, '.mat'])
        save('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'time_res', 'timestamp')
end


%% fit polarization


% % phase change
%     new_scan_list   = mean(scan_list_rbv, 2)-mean(scan_list_set);
%     new_scan_list   = scan_list_set-mean(scan_list_set);
%     
%    
% 
%     % x fit
%     pos_CoM_x_mean  =  scale_x * ( mean(pos_CoM_x, 2) - 0*mean(pos_CoM_x(:))) ;  % pixel
%     pos_CoM_x_std   =  scale_x * std(pos_CoM_x, 1, 2); 
%  
%     % y fit
%     pos_CoM_y_mean  =  scale_y * ( mean(pos_CoM_y, 2) - 0*mean(pos_CoM_y(:))) ;  % pixel
%     pos_CoM_y_std   =  scale_y * std(pos_CoM_y, 1, 2); 
%     
%     [par_xy, yFit_xy, parstd_xy] = util_polyFit(pos_CoM_x_mean, pos_CoM_y_mean, 1);
%     
%     polarization_TDS     = rad2deg(atan(par_xy(1))); % in deg
%     polarization_TDS_err = rad2deg(parstd_xy(1)/(par_xy(1)^2+1));
%     
%     
%     % r list ----> needs update !!!!
% %     pos_CoM_r       =  sqrt( (pos_CoM_x - pos_CoM_x(10) ).^2 + ( pos_CoM_y - pos_CoM_y(10) ).^2); 
%     pos_CoM_r       =  sqrt( (pos_CoM_x-min(pos_CoM_x(:)) ).^2 + ( pos_CoM_y-min(pos_CoM_y(:)) ).^2);
%     pos_CoM_r_mean  =  mean(pos_CoM_r, 2);
%     pos_CoM_r_std   =  std(pos_CoM_r, 1, 2);  
%     [par_r, yFit_r, parstd_r] = util_polyFit(new_scan_list, pos_CoM_r_mean, 1, pos_CoM_r_std);
%     
%     timecal_fspixel = 1e15 * deg2rad(1/par_r(1))/(2*pi*frequency_XTDS);
%     timecal_fspixel_err = 1e15 * deg2rad(parstd_r(1)/par_r(1)^2)/(2*pi*frequency_XTDS);
%     
%     timecal_streak  = 1/(3e8/(1e-3*scale_x) * deg2rad(1/par_r(1))/(2*pi*frequency_XTDS));
%   %  timecal_streak_err  = 3e8/(1e-3*scale_x) * deg2rad(parstd_r(1)/par_r(1)^2)/(2*pi*frequency_XTDS);
%     
%     
% % final plot of polarization -- expert mode only
% 
% figure(5)
% set(gcf, 'OuterPosition', [100, 100, 1000, 500])
% 
%     subplot(1,2,1)
% 
%     p1 = plot(scale_x * pos_CoM_x(:),scale_y * pos_CoM_y(:), 'b.');
%     hold on
%         p2 = plot(pos_CoM_x_mean, yFit_xy, '-r');
%     hold off
%     grid on
%     axis equal
%     title(['time calibration ', name_cam], 'Interpreter', 'none', 'FontSize', fontSize)
%     xlabel(['beam position <x> (mm)'], 'FontSize', fontSize)
%     ylabel(['beam position <y> (mm)'], 'FontSize', fontSize)
%     legend([p1(1), p2], ...
%     {'data', ['fit pol. angle: ', num2str(polarization_TDS, '%5.2f'), ' +/- ', num2str(polarization_TDS_err, '%5.2f'), ' deg']}, ...
%     'Location', 'NorthWest', 'FontSize', fontSize-1)
%     set(gca, 'FontSize', fontSize) 
%     
%     subplot(1,2,2)
%     p1 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*pos_CoM_r, 'b.');
%     hold on
%         p2 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*yFit_r, '-r');
%     hold off
%     grid on
%     title(timestamp, 'Interpreter', 'none', 'FontSize', fontSize)
%     xlabel([name_actuator, ' (ps)'], 'FontSize', fontSize)
%     ylabel('beam position sqrt(<x>^2+<y>^2) (mm)', 'FontSize', fontSize) 
%     legend([p1(1), p2], ...
%         {'data', ['fit time cal: ', num2str(timecal_fspixel, '%5.2f'), ' +/- ', num2str(timecal_fspixel_err, '%5.2f'), ' fs/pixel', 10 ,...
%                            '(fit streak: ', num2str(timecal_streak, '%5.1f'), ')']}, ...
%         'Location', 'NorthWest', 'FontSize', fontSize-1)
%     set(gca, 'FontSize', fontSize) 
%     
%      % add elog printing button
%      uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX time calibration', comment));
% 
% 
% 
% 
% 
% %% save 
% 
% % if (flag_save)
% %         save(['/home/ttflinac/user/swesch/PolariX/time_calibration/', name_script(5:end), '_', timestamp, '.mat'])
% % end
% 
% if (flag_save)
%         save(['/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/', name_script(5:end), '_', timestamp, '.mat'])
%         save('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_polarization_9FL2XTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'timestamp')
% end


end
