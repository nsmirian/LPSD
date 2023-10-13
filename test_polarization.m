%% path

    path_tools          = '/home/ttflinac/user/swesch/Tools';
    addpath(genpath(path_tools));
    col                 = gen_ColorDefinitions;
    load('myColorMap.mat') % cmap

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
    time_res = (sigma_x * 1e-3 / 3e8 ) / timecal_streak *1e15 ; % in m
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
                           '(streak ', num2str(timecal_streak, '%5.1f'), ')', 10, ...
                           'Time res: ' num2str((time_res), '%5.1f'), ' fs']}, ...
        'Location', 'NorthWest', 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize) 
    
     % add elog printing button
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FLF PolariX time calibration', comment));


     % reset camera
%     tmp      = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 0);
%     tmp      = doocswrite([addr_cam, 'EXPOSURETIMEABS'], 1000);


%% save 

if (flag_save)
        save(['/home/ttflinac/user/swesch/PolariX/FL2/time_calibration/', name_script(5:end), '_', timestamp, '.mat'])
        save('/home/ttflinac/user/swesch/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'time_res', 'timestamp')
end


%% fit polarization


% phase change
    new_scan_list   = mean(scan_list_rbv, 2)-mean(scan_list_set);
    new_scan_list   = scan_list_set-mean(scan_list_set);
    
   

    % x fit
    pos_CoM_x_mean  =  scale_x * ( mean(pos_CoM_x, 2) - 0*mean(pos_CoM_x(:))) ;  % pixel
    pos_CoM_x_std   =  scale_x * std(pos_CoM_x, 1, 2); 
 
    % y fit
    pos_CoM_y_mean  =  scale_y * ( mean(pos_CoM_y, 2) - 0*mean(pos_CoM_y(:))) ;  % pixel
    pos_CoM_y_std   =  scale_y * std(pos_CoM_y, 1, 2); 
    
    [par_xy, yFit_xy, parstd_xy] = util_polyFit(pos_CoM_x_mean, pos_CoM_y_mean, 1);
    
    polarization_TDS     = rad2deg(atan(par_xy(1))); % in deg
    polarization_TDS_err = rad2deg(parstd_xy(1)/(par_xy(1)^2+1));
    
    
    % r list ----> needs update !!!!
%     pos_CoM_r       =  sqrt( (pos_CoM_x - pos_CoM_x(10) ).^2 + ( pos_CoM_y - pos_CoM_y(10) ).^2); 
    pos_CoM_r       =  sqrt( (pos_CoM_x-min(pos_CoM_x(:)) ).^2 + ( pos_CoM_y-min(pos_CoM_y(:)) ).^2);
    pos_CoM_r_mean  =  mean(pos_CoM_r, 2);
    pos_CoM_r_std   =  std(pos_CoM_r, 1, 2);  
    [par_r, yFit_r, parstd_r] = util_polyFit(new_scan_list, pos_CoM_r_mean, 1, pos_CoM_r_std);
    
    timecal_fspixel = 1e15 * deg2rad(1/par_r(1))/(2*pi*frequency_XTDS);
    timecal_fspixel_err = 1e15 * deg2rad(parstd_r(1)/par_r(1)^2)/(2*pi*frequency_XTDS);
    
    timecal_streak  = 1/(3e8/(1e-3*scale_x) * deg2rad(1/par_r(1))/(2*pi*frequency_XTDS));
  %  timecal_streak_err  = 3e8/(1e-3*scale_x) * deg2rad(parstd_r(1)/par_r(1)^2)/(2*pi*frequency_XTDS);
    
    
%% final plot

figure(5)
set(gcf, 'OuterPosition', [100, 100, 1000, 500])

    subplot(1,2,1)

    p1 = plot(scale_x * pos_CoM_x(:),scale_y * pos_CoM_y(:), 'b.');
    hold on
        p2 = plot(pos_CoM_x_mean, yFit_xy, '-r');
    hold off
    grid on
    axis equal
    title(['time calibration ', name_cam], 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel(['beam position <x> (mm)'], 'FontSize', fontSize)
    ylabel(['beam position <y> (mm)'], 'FontSize', fontSize)
    legend([p1(1), p2], ...
    {'data', ['fit pol. angle: ', num2str(polarization_TDS, '%5.2f'), ' +/- ', num2str(polarization_TDS_err, '%5.2f'), ' deg']}, ...
    'Location', 'NorthWest', 'FontSize', fontSize-1)
    set(gca, 'FontSize', fontSize) 
    
    subplot(1,2,2)
    p1 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*pos_CoM_r, 'b.');
    hold on
        p2 = plot(1e12*deg2rad(new_scan_list)/(2*pi*frequency_XTDS), scale_x*yFit_r, '-r');
    hold off
    grid on
    title(timestamp, 'Interpreter', 'none', 'FontSize', fontSize)
    xlabel([name_actuator, ' (ps)'], 'FontSize', fontSize)
    ylabel('beam position sqrt(<x>^2+<y>^2) (mm)', 'FontSize', fontSize) 
    legend([p1(1), p2], ...
        {'data', ['fit time cal: ', num2str(timecal_fspixel, '%5.2f'), ' +/- ', num2str(timecal_fspixel_err, '%5.2f'), ' fs/pixel', 10 ,...
                           '(fit streak: ', num2str(timecal_streak, '%5.1f'), ')']}, ...
        'Location', 'NorthWest', 'FontSize', fontSize-1)
    set(gca, 'FontSize', fontSize) 
    
     % add elog printing button
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FLF PolariX time calibration', comment));