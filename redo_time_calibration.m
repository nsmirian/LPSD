%% redo calibration
load('/home/ttflinac/user/swesch/PolariX/FL2/time_calibration/timecalibration_OTR9FL2XTDS_reso_2021-04-05T164835.mat')

%% path

    path_tools          = '/home/ttflinac/user/swesch/Tools';
    addpath(genpath(path_tools));
    col                 = gen_ColorDefinitions;
    load('myColorMap.mat') % cmap

%% pre allocating    

    
    % remove scan points
    num_actuator = 5;
    start=1;
    scan_list_set = scan_list_set(start:start+num_actuator-1);
    scan_list_rbv = scan_list_rbv(start:start+num_actuator-1);
    
            % center of mass
    pos_CoM_x       = zeros(num_actuator, num_sig);
    pos_CoM_y       = zeros(num_actuator, num_sig);
%% analyze images
figure;
for ii = 1:num_actuator % scan points
    
       
     
    %%%%%%%%%%  analyze data %%%%%%%%%%

    for jj = 1:num_sig
        

        % remove bg x
%         tmp                     = squeeze(raw_spec_x(ii,jj,:));
%         corr_spec_x(ii,jj,:)    = tmp' - bgr_spec_x_mean;
 
        % filter (mean), fit (asymmetric gauss) and position (from fit)
        tmp                     = medfilt1(corr_spec_x(ii,jj,:), 3);
        clf
        plot(squeeze(tmp))
        hold on
        [tmp, yFit]             = util_gaussFit(1:length(tmp), tmp, 1, 1);
        plot(squeeze(yFit))
        plot(round(tmp(2)), 5, 'x', 'MarkerSize', 20, 'LineWidth', 4)
        pos_CoM_x(ii,jj)        = tmp(2);
        
        
        
        % remove bg y
%         tmp                     = squeeze(raw_spec_y(ii,jj,:));
%         corr_spec_y(ii,jj,:)    = tmp' - bgr_spec_y_mean;
        
        % filter (mean), fit (asymmetric gauss) and position (from fit)
        tmp                     = medfilt1(corr_spec_y(ii,jj,:), 3);
        [tmp, tmp2]             = util_gaussFit(1:length(tmp), tmp, 1, 1);
        pos_CoM_y(ii,jj)        = tmp(2);        
        
        
        
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
    
    time_res = (sigma_x * 1e-3 / 3e8 ) / timecal_streak *1e15 ; % in m
    
    %% final plot
figure(5)
clf
plot(pos_CoM_x.')
figure(4)
clf
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


