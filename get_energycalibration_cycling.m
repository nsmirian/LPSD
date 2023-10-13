
%% script for PolariX energy calibration for screen OTR9FL2XTDS
%
% swesch - 8th Sep 2019 - v1: first script 
% 
%                         v2: 

% clear all
% close all

function get_energycalibration_cycling(start_actuator, end_actuator, block_laser, num_bgr, num_sig, num_actuator, comment)    
%% config

%     comment         = 'FL2 PolariX energy calibration';

    flag_save       = 1;
    flag_act        = 1;
    flag_pause      = 1; % make a pause between reach goal current and take measurement

    % D4FL2XTDS current range
%     start_actuator  = 117; % in A 
%     end_actuator    = 121; % in A 

%     num_actuator    = 5;    % set points
%     num_bgr         = 10;    % number of single background images per set point
%     num_sig         = 10;    % number of single measurement per set point

% 	block_laser = 0; % 1 = use injector laser; 0 = use trigger delay

    fontSize        = 14;


%% 

    % addresses
    name_script             = 'get_energycalibration_Cycling';
    name_cam                = 'OTR9FL2XTDS';
    addr_cam                = ['FLASH.DIAG/CAMERA/', name_cam, '/'];
    delay = doocsread([addr_cam, 'TRIGGERDELAYABS']);
    
    pixelformat=doocswrite([addr_cam, 'FORMAT.OUT'], 0); %until IAS is used: set pixelformat out to mono8

    % address magnet directly
    name_magnet             = 'D4FL2XTDS';
    addr_actuator_set       = ['FLASH.MAGNETS/MAGNET.ML/', name_magnet, '/CURRENT.SP'];
    addr_actuator_rbv       = ['FLASH.MAGNETS/MAGNET.ML/', name_magnet, '/CURRENT.RBV'];

    % time
    timestamp               = datestr(clock, 'yyyy-mm-ddTHHMMSS');


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

    % current rbv
    actuator_rbv    = zeros([1, num_actuator]);

    % charge
    charge_7FLFMAFF = zeros(num_actuator, num_sig);
    charge_7FLFDUMP = zeros(num_actuator, num_sig);

%% prepare block laser for FLASH2

    % which laser
    name_laser          = getfield(doocsread('FLASH.DIAG/TIMER/FLASHCPUTIME1.0/LASER_SELECT.2'), 'data');
%     addr_laser_block    = ['FLASH.DIAG/LASER.CONTROL/LASER', num2str(name_laser), '/BLOCK_LASER'];
    % addr_laser_block    = ['FLASH.DIAG/BEAMLINES/FLASH/LASER/BLOCK_LASER.FLASH3'];
    
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

%% Cycling befor measerment

ddd_write = doocswrite(addr_laser_block, 1);  % block laser first , 
ddd_write = doocswrite(addr_actuator_set,start_actuator);
% now cycling 
addr_actuator_Cycling       = ['FLASH.MAGNETS/MAGNET.ML/', name_magnet, 'PS.CYCLE'];

ddd_write = doocswrite(addr_actuator_Cycling ,1);

% wait till cycling finish

t=doocsread('TTF2.MAGNETS/QUAD/Q22FL2EXTR/PS.CLEAN'); 
cleanV=t.data;
while cleanV                                        %% we should aask about this part 
    display_message(msgLoc, [' Dipole is cycling, please wait']); 
    pause(60)
end
display_message(msgLoc, [' Dipole cycling finished ']); 

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
        title(['background image SCR', name_cam], 'FontSize', fontSize)
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

for ii = 1:num_actuator % scan points
    
       
    %%%%%%%%%% set dipole current %%%%%%%%%%     %%%%%
    if (flag_act)
                        
        % set value 
        ddd_write = doocswrite(addr_actuator_set, scan_list_set(ii));
        display_message(msgLoc, [name_script, '(): Set actuator ', name_magnet, ' to ' , num2str(scan_list_set(ii), '%5.3f')]);
        
        % wait for set = rbv
        display_message(msgLoc, [name_script, '(): Actuator ', name_magnet, ' set.']);
        
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

%%
% restore actuator reference value 
if (flag_act)
    tmp = mean([start_actuator, end_actuator]);
    ddd_write = doocswrite(addr_actuator_set, tmp);
    display_message(msgLoc, [name_script, '(): Set actuator ', name_magnet, ' to ' , num2str(tmp, '%5.3f')]);
end

% Cycling befor measerment

ddd_write = doocswrite(addr_laser_block, 1);  % block laser first , 

% now cycling 
addr_actuator_Cycling       = ['FLASH.MAGNETS/MAGNET.ML/', name_magnet, 'PS.CYCLE'];

ddd_write = doocswrite(addr_actuator_Cycling ,1);

% wait till cycling finish
t=doocsread('TTF2.MAGNETS/QUAD/Q22FL2EXTR/PS.CLEAN'); 
cleanV=t.data;
while CleanV                                       %   we should aask about this part 
    display_message(msgLoc, [' Dipole is cycling, please wait']); 
    pause(60)
end
display_message(msgLoc, [' Dipolecycling finished ']); 

%% clear stuff
clear tmp tmp2 ddd_write ddd_read 

%% fit 

    % compute rel. current change
    rel_scan_list   = mean(scan_list_rbv, 2)/mean(scan_list_set)-1;

    % x fit
    pos_CoM_x_mean  =  mean(pos_CoM_x, 2);  % pixel
    pos_CoM_x_std   =  std(pos_CoM_x, 1, 2); 
    [par_x, yFit_x, parstd_x] = util_polyFit(rel_scan_list, pos_CoM_x_mean, 1, pos_CoM_x_std);

    % y fit
    pos_CoM_y_mean  =  mean(pos_CoM_y, 2);  % pixel
    pos_CoM_y_std   =  std(pos_CoM_y, 1, 2); 
    [par_y, yFit_y, parstd_y] = util_polyFit(rel_scan_list, pos_CoM_y_mean, 1, pos_CoM_y_std);
    
    
    ergcal          = 1e3/par_y(1); % permille per pixel
    ergcal_err      = 1e3*parstd_y(1)/par_y(1)^2;
 
    
    
%% final plot

figure(4)
set(gcf, 'OuterPosition', [100, 100, 1000, 500])

    subplot(1,2,1)

    p1 = plot(rel_scan_list, pos_CoM_y, 'b.');
    hold on
        p2 = plot(rel_scan_list, yFit_y, '-r');
    hold off
    grid on
    title(['energy calibration ', name_cam], 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel(['relative ', name_magnet, ' current'], 'FontSize', fontSize)
    ylabel([name_cam, ' beam position <x> (pixel)'], 'FontSize', fontSize)
    legend([p1(1), p2], ...
            {'data', ['fit: ', num2str(1e3/par_y(1), '%5.3f'), ' +/- ', num2str(1e3*parstd_y(1)/par_y(1)^2, '%5.3f'), ' 1e-3/pixel', ...
            num2str(1e3/par_y(1)/scale_x, '%5.3f'), ' +/- ', num2str(1e3*parstd_y(1)/par_y(1)^2/scale_x, '%5.3f'), ' 1e-3/mm']}, ...
            'Location', 'NorthWest', 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize) 
    
    subplot(1,2,2)
    p1 = plot(1e3*rel_scan_list, scale_y*pos_CoM_y, 'b.');
    hold on
        p2 = plot(1e3*rel_scan_list, scale_y*yFit_y, '-r');
    hold off
    grid on
    title(timestamp, 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel(['relative ', name_magnet, ' current / 10^{-3}'], 'FontSize', fontSize)
   % ylabel('beam position <x> (mm)', 'FontSize', fontSize) 
    legend([p1(1), p2], ...
        {'data', ['fit: ', num2str(1e-3*scale_y*par_y(1), '%5.4f'), ' +/- ', num2str(1e-3*scale_y*parstd_y(1), '%5.4f'), ' m']}, ...
        'Location', 'NorthWest', 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize) 
    
     % add elog printing button
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FLF PolariX energy calibration', comment));





%% save 

if (flag_save)
        save(['/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/', name_script(5:end), '_', timestamp, '.mat'])
        save('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err', 'timestamp')
end


end