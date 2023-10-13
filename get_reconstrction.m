%% take image from abr. beam camera in FL2 
% 
%  make sure that: (1) camera is not saturated
%                  (2) 
%  Author: Najmeh Mirian 
% date=30 Agu 2023
% Flash control room
% reconstruction 



% 
% clear all
% close all
 
%    
function get_calibratedimage_OTR9FL2XTDS(block_laser, num_bgr, num_sig,  comment)   
% comment         = 'FL2 PolariX measurement';
flag_save       = 1;
flag_act        = 1;

%
% num_bgr         = 10;  % number of single background images per set point
% num_sig         = 10;  % number of single measurement per set point

%
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');
timestamp1=timestamp ;
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

    %display_message(' - setup done');

%% take background 
    
% block laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 1);
        %display_message(msgLoc, [' - block LASER', name_laser]);  
        pause(2)         
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 10000);
   % display_message( [' - changed camera delay']);  
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
            %display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
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


%display_message(' - background data taken');

%% take data


% unblock laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 0); % unblock laser 
        %display_message(msgLoc, [name_script, '(): Unblock LASER', name_laser]);
        pause(0.5)       
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], delay.data);
 %   display_message( [' - changed camera delay back']);  
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
        while img_sig(:,:,jj) == img_sig(:,:,jj-1)                              % change it to time stamp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ddd_read        = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            img_sig(:,:,jj) = ddd_read.data.val_val;   
       %     display_message( [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(1/rep_rate)
        end
    end      
    
    % read charge in nC
    ddd_read = doocsread('FLASH.DIAG/TOROID/7FL2XTDS/CHARGE.FLASH2');
    charge_7FL2XTDS(jj) = ddd_read.data;

    pause(2/rep_rate)

end
    
%display_message(msgLoc, ' - data taken');

% compute mean()
charge_7FL2XTDS_mean = mean(charge_7FL2XTDS); 

% clear stuff
clear tmp ddd_write ddd_read



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
y=1:length_y;
%%
for jj = 1:num_sig
    
   tmp_img             = hlc_clean_image( squeeze(img_sig(:,:,jj)) - img_bgr_mean );
    %tmp_img             = get_ROI(squeeze(img_sig(:,:,jj)) - img_bgr_mean);
    img_filt(:,:,jj)    = medfilt2(tmp_img);
    
    % time
    tmp_profile         = mean(img_filt(:,:,jj));
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
cent_energy_cr=zeros(num_sig, length_x);
slice_enrgy_spread=zeros(num_sig, length_x);
%%
% tmp2_img = squeeze(transpose(img_filt(:,:,end)));
% % jj=1
% 
%     for n=min(time_pos_good):max(time_pos_good)
%               
%                 mu0=  find(tmp2_img(n,erg_pos_good)==max(tmp2_img(n,erg_pos_good)));
%                 [sigma, mu] = gaussfit(erg_pos_good, smoothdata(tmp2_img(n,erg_pos_good)) , 0, erg_pos_good(mu0(1)));
%                 if isnan(sigma)
%                 else
% 
%                 slice_enrgy_spread_ave( n)=sigma*calib_y;
%                 cent_energy_cr_ave( n)=mu*calib_y;
%                 end
%     end
% 
%     figure; plot(slice_enrgy_spread_ave)
%%
for jj=1:num_sig
    for n=min(time_pos_good):max(time_pos_good)
              
                mu0=  find( img_filt(erg_pos_good,n,jj)==max(img_filt(erg_pos_good,n,jj)));
                [sigma, mu] = gaussfit( erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) , 0, erg_pos_good(mu0(1)));
                if isnan(sigma)
                else

                slice_enrgy_spread(jj, n)=sigma*calib_y;
                cent_energy_cr(jj, n)=mu*calib_y;
                end
                 
    end
end


%% jetting of stabilitty correction 
for jj=1: num_sig
        mass_center(jj)=find(x_profile(jj,:)==max(x_profile(jj,:)));
end
%%
figure(1044);hold on

sumes=zeros( length_x);
sumEcr=zeros( length_x);
sumCu=zeros( length_x);

for jj=1:num_sig
    subplot(3,3,1);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),slice_enrgy_spread(jj,time_pos_good-mass_center(1)+mass_center(jj)))
    sumes(time_pos_good)=sumes(time_pos_good)+slice_enrgy_spread(jj,time_pos_good-mass_center(1)+mass_center(jj));
    subplot(3,3,2);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),cent_energy_cr(jj,time_pos_good-mass_center(1)+mass_center(jj)))
    sumEcr(time_pos_good)=sumEcr(time_pos_good)+cent_energy_cr(jj,time_pos_good-mass_center(1)+mass_center(jj));
    subplot(3,3,3);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),1e-3*charge_7FL2XTDS(jj)*1e-9*1e15*x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj)))
    sumCu(time_pos_good)=sumCu(time_pos_good)+x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj));
end
average_slice_energy_spread=sumes/num_sig;
aver_central_energy=sumEcr/num_sig;
aver_current=sumCu/num_sig*1e-3*charge_7FL2XTDS_mean*1e-9;

subplot(3,3,1)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), average_slice_energy_spread(time_pos_good),'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
xlabel('time t /fs', 'FontSize', fontSize)        
        ylabel('slice energy Sepread (Mev)', 'FontSize', fontSize)
subplot(3,3,2)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_central_energy(time_pos_good),'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
xlabel('time t /fs', 'FontSize', fontSize)        
ylabel('\Delta energy cneter (Mev)', 'FontSize', fontSize)
subplot(3,3,3)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_current(time_pos_good),'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
xlabel('time t /fs', 'FontSize', fontSize)        
ylabel('peak current /kA', 'FontSize', fontSize)
%%

% disp(' - ')
figure(1)
set(gcf, 'Position', [400, 1, 700, 700])

    subplot(2,4,1)
     imagesc(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            erg_axis(erg_pos_good), flipud(tmp2_img(time_pos_good,erg_pos_good)') )
        grid on
        title(['single image ', 'with SASE'], 'FontSize', fontSize, 'Interpreter', 'none')
        xlabel('time t /fs)', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}' , 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize,'YDir', 'normal')
        lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
        
    subplot(2,4,2)       
        p1 = plot(y_profile(end,erg_pos_good), -erg_axis(erg_pos_good) );
        grid on
        title(timestamp1, 'FontSize', fontSize)
        xlabel('norm. erg density', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}', 'FontSize', fontSize)
        legend([p1(1)], {['rms: ' num2str(y_var(end), '%5.2f'), ' \cdot 10^{-3}', 10, ...
                'FWHM: ', num2str(y_fwhm(end), '%5.2f'), ' \cdot 10^{-3}']}, ...
                'FontSize', fontSize-2)        
        set(gca, 'FontSize', fontSize, 'YDir', 'normal');%, 'Ylim', -lim_y)        
        
    subplot(2,4,3)
  
    yyaxis left
        p1 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            1e-3*charge_7FL2XTDS_mean*1e-9*1e15 * x_profile(end,time_pos_good));
        grid on
        title('current and energy spread', 'Interpreter', 'none')
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
        yyaxis right
        
           p2 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), average_slice_energy_spread(time_pos_good))
             
           ylabel('energy spread (Mev)')
        legend([p2(1),p1(1)], { ['energy spread'],['current   ','rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
                        'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
                        'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
                        'FontSize', fontSize-2)
                

        set(gca, 'FontSize', fontSize, 'Xlim', lim_x)           
                %'FWHM: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
        %         

    subplot(2,4,4)
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
    % filename=['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp1, '.mat']
   %  text=[comment '\r\n', ...
        %     'data saved in:', filename]

     %comment='data saved at '
   %  uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX LPSD@OTR9FL2XTDS', text));
   %% save 

if (flag_save)
        save(['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp1,'_', 'with_SASE' , '.mat'], ...
             'comment', 'timestamp', 'name_cam', 'img_bgr', 'img_sig', 'charge_7FL2XTDS', ...
             'timecal_fspixel', 'timecal_fspixel_err', 'timeres_calib', 'ergcal', 'ergcal_err', 'phase_TDS', 'tmp2_img', ...
              'x_com', 'x_var', 'x_fwhm', 'x_axis', 'x_profile',...
              'y_com', 'y_var', 'y_fwhm', 'y_axis', 'y_profile');
        %display_message([' - data saved']);
end


     %%   Second data   without SASE  kill sase by kick 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear img_bgr img_sig img_bgr_mean
launchfb_add='FLASH.FEEDBACK/ORBIT.FLASH2/ORBITFEEDBACK/ACTIVATE_FB';
ddd_write = doocswrite(launchfb_add, 0);


lunch_H3FL2SEED4=getfield(doocsread('FLASH.MAGNETS/MAGNET.ML/H3FL2SEED4/CURRENT.SP'), 'data');
new_lunch_H3FL2SEED4=doocswrite('FLASH.MAGNETS/MAGNET.ML/H3FL2SEED4/CURRENT.SP',lunch_H3FL2SEED4+0.3 );   % we killed SASE

timestamp2       = datestr(clock, 'yyyy-mm-ddTHHMMSS');
%% take background 
    
% block laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 1);
   %     display_message(msgLoc, [' - block LASER', name_laser]);  
        pause(2)         
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], 10000);
   % display_message( [' - changed camera delay']);  
    pause(2)
end

% take bgr
for jj = 1:num_bgr 
    
    % read
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_bgr_no_sase(:,:,jj)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_bgr_no_sase(:,:,jj) == img_bgr_no_sase(:,:,jj-1)
            ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            img_img_bgr_no_sase(:,:,jj)       = ddd_read.data.val_val;
            %display_message(msgLoc, [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(10/rep_rate)
        end
    end
    
    pause(1/rep_rate)
    
end
    
% compute mean()
if num_bgr == 1
    img_bgr_mean_no_sase = squeeze(img_bgr_no_sase);
else
    img_bgr_mean_no_sase = squeeze(mean(img_bgr_no_sase, 3));
end


%display_message(' - background data taken');

%% take data


% unblock laser
if block_laser
    if (flag_act)
        ddd_write = doocswrite(addr_laser_block, 0); % unblock laser 
        %display_message(msgLoc, [name_script, '(): Unblock LASER', name_laser]);
        pause(0.5)       
    end
else
    ddd_write = doocswrite([addr_cam, 'TRIGGERDELAYABS'], delay.data);
    %display_message( [' - changed camera delay back']);  
    pause(2)
end


   
% loop
charge_7FL2XTDS_no_sase     = zeros([1, num_sig]);
phase_TDS_no_sase           = getfield(doocsread('FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE'),'data');
for jj = 1:num_sig

    %%% read spectrum x
    ddd_read            = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
    img_sig_no_sase(:,:,jj)     = ddd_read.data.val_val;
    
    % compared it with last measurement
    if jj > 1
        while img_sig_no_sase(:,:,jj) == img_sig_no_sase(:,:,jj-1)                              % change it to time stamp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ddd_read        = doocsread([addr_cam, 'IMAGE_EXT_ZMQ']);
            img_sig_no_sase(:,:,jj) = ddd_read.data.val_val;   
            %display_message( [' - (', num2str(jj), ') same data ... wait ...']);  
            pause(1/rep_rate)
        end
    end      
    
    % read charge in nC
    ddd_read = doocsread('FLASH.DIAG/TOROID/7FL2XTDS/CHARGE.FLASH2');
    charge_7FL2XTDS_no_sase(jj) = ddd_read.data;

    pause(2/rep_rate)

end
    
%display_message(msgLoc, ' - data taken');

% compute mean()
charge_7FL2XTDS_mean_no_sase = mean(charge_7FL2XTDS_no_sase); 

% clear stuff
clear tmp ddd_write ddd_read img_filt tmp_img tmp2_img mass_center





%lunch_H3FL2SEED4=getfield(doocsread('FLASH.MAGNETS/MAGNET.ML/H3FL2SEED4/CURRENT.SP'), 'data');
new_lunch_H3FL2SEED4=doocswrite('FLASH.MAGNETS/MAGNET.ML/H3FL2SEED4/CURRENT.SP',lunch_H3FL2SEED4 );   % we killed SASE
ddd_write = doocswrite(launchfb_add, 1);

%% 

% analysis
img_filt_no_sase    = img_sig_no_sase;

num_sig     = size(img_sig_no_sase, 3);
length_y    = size(img_sig_no_sase, 1);
length_x    = size(img_sig_no_sase, 2);
img_bgr_mean = squeeze(mean(img_bgr, 3));

x_com       = zeros([1, num_sig]);              y_com       = zeros([1, num_sig]);
x_var       = zeros([1, num_sig]);              y_var       = zeros([1, num_sig]);
x_fwhm      = zeros([1, num_sig]);              y_fwhm      = zeros([1, num_sig]);
x_profile   = zeros([num_sig, length_x]);       y_profile   = zeros([num_sig, length_y]);
y=1:length_y;
%%
for jj = 1:num_sig
    
   tmp_img             = hlc_clean_image( squeeze(img_sig_no_sase(:,:,jj)) );%- 0*img_bgr_mean_no_sase );
  %  tmp_img             = get_ROI(squeeze(img_sig(:,:,jj)) - img_bgr_mean);
    img_filt_no_sase(:,:,jj)    = tmp_img;
    
    % time
    tmp_profile         = mean(tmp_img);
%     tmp_profile(2010)=tmp_profile(2009); % remove hot pixel
%     tmp_profile         = medfilt1(tmp_profile, 1);
    [x_com(jj), x_var(jj), x_fwhm(jj), x_axis, x_profile(jj,:)]  = get_profile_stats(tmp_profile, calib_x);
    
    % enrg
    tmp_profile         = mean(tmp_img, 2);
%     tmp_profile         = medfilt1(tmp_profile, 1);
    [y_com(jj), y_var(jj), y_fwhm(jj), y_axis, y_profile(jj,:)]  = get_profile_stats(tmp_profile', calib_y);

              
   
end
%% slice analysis


fontSize        = 14;

erg_pos_good =  find(y_profile(end,:)>0); erg_pos_good = min(erg_pos_good):max(erg_pos_good);
time_pos_good = find(x_profile(end,:)>0); time_pos_good = min(time_pos_good):max(time_pos_good);
cent_energy_cr=zeros(num_sig, length_x);
slice_enrgy_spread=zeros(num_sig, length_x);
%%
% tmp2_img = squeeze(transpose(img_filt_no_sase(:,:,end)));
% % jj=1
% 
%     for n=min(time_pos_good):max(time_pos_good)
%               
%                 mu0=  find(tmp2_img(erg_pos_good,n)==max(tmp2_img(erg_pos_good,n)));
%                 [sigma, mu] = gaussfit( erg_pos_good, smoothdata(tmp2_img(erg_pos_good,n)) , 0, erg_pos_good(mu0(1)));
%                 if isnan(sigma)
%                 else
% 
%                 slice_enrgy_spread_ave( n)=sigma*calib_y;
%                 cent_energy_cr_ave( n)=mu*calib_y;
%                 end
%     end
% 
%     figure; plot(slice_enrgy_spread_ave)
%%
for jj=1:num_sig
    for n=min(time_pos_good):max(time_pos_good)
              
                mu0=  find( img_filt(erg_pos_good,n,jj)==max(img_filt(erg_pos_good,n,jj)));
                [sigma, mu] = gaussfit( erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) , 0, erg_pos_good(mu0(1)));
                if isnan(sigma)
                else

                slice_enrgy_spread_nosase(jj, n)=sigma*calib_y;
                cent_energy_cr_nosase(jj, n)=mu*calib_y;
                end
                 
    end
end



%% jetting of stabilitty correction 
for jj=1: num_sig
mass_center(jj)=find(x_profile(jj,:)==max(x_profile(jj,:)));
end
%%
figure(1044);hold on

sumes_nosase=zeros( length_x);
sumEcr_nosase=zeros( length_x);
sumCu_nosase=zeros( length_x);

for jj=1:num_sig
    subplot(3,3,4);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),slice_enrgy_spread_nosase(jj,time_pos_good-mass_center(jj)+mass_center(jj)))
    sumes_nosase(time_pos_good)=sumes_nosase(time_pos_good)+slice_enrgy_spread_nosase(jj,time_pos_good-mass_center(1)+mass_center(jj));
    subplot(3,3,5);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),cent_energy_cr_nosase(jj,time_pos_good-mass_center(1)+mass_center(jj)))
    sumEcr_nosase(time_pos_good)=sumEcr_nosase(time_pos_good)+cent_energy_cr_nosase(jj,time_pos_good-mass_center(1)+mass_center(jj));
    subplot(3,3,6);hold on
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),1e-3*charge_7FL2XTDS_no_sase(jj)*1e-9*1e15*x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj)))
    sumCu_nosase(time_pos_good)=sumCu_nosase(time_pos_good)+x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj));
end
average_slice_energy_spread_nosase=sumes_nosase /num_sig;
aver_central_energy_nosase        =sumEcr_nosase/num_sig;
aver_current_nosase               =sumCu_nosase/num_sig*1e-3*charge_7FL2XTDS_mean_no_sase*1e-9*1e15;
subplot(3,3,4)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), average_slice_energy_spread_nosase(time_pos_good),'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
xlabel('time t /fs', 'FontSize', fontSize)        
        ylabel('slice energy Sepread (Mev)', 'FontSize', fontSize)
subplot(3,3,5)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_central_energy_nosase(time_pos_good) ,'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
xlabel('time t /fs', 'FontSize', fontSize)        
ylabel('\Delta energy cneter (Mev)', 'FontSize', fontSize)
subplot(3,3,6)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_current_nosase(time_pos_good),'--k', 'LineWidth',3 )

set(gca, 'FontSize', fontSize,'YDir', 'normal')
lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
xlabel('time t /fs', 'FontSize', fontSize)        
ylabel('peak current /kA', 'FontSize', fontSize)


%% plot 1




tmp2_img = squeeze(transpose(img_filt_no_sase(:,:,end)));

% disp(' - ')
figure(1)
set(gcf, 'Position', [400, 1, 700, 700])

    subplot(2,4,5)
        imagesc(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            erg_axis(erg_pos_good), flipud(tmp2_img(time_pos_good,erg_pos_good)') )
        grid on
        title(['single image No SASE ', name_cam], 'FontSize', fontSize, 'Interpreter', 'none')
        xlabel('time t /fs)', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}' , 'FontSize', fontSize)
        set(gca, 'FontSize', fontSize,'YDir', 'normal')
        lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
        
    subplot(2,4,6)       
        p1 = plot(y_profile(end,erg_pos_good), -erg_axis(erg_pos_good) );
        grid on
        title('energy profile', 'FontSize', fontSize)
        xlabel('norm. erg density', 'FontSize', fontSize)        
        ylabel('rel. energy deviation \delta / 10^{-3}', 'FontSize', fontSize)
        legend([p1(1)], {['rms: ' num2str(y_var(end), '%5.2f'), ' \cdot 10^{-3}', 10, ...
                'FWHM: ', num2str(y_fwhm(end), '%5.2f'), ' \cdot 10^{-3}']}, ...
                'FontSize', fontSize-2)        
        set(gca, 'FontSize', fontSize, 'YDir', 'normal');%, 'Ylim', -lim_y)        
        
    subplot(2,4,7)
    
    yyaxis left
        p1 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
            1e-3*charge_7FL2XTDS_mean_no_sase*1e-9*1e15 * x_profile(end,time_pos_good));
        grid on
        title('current and energy spread', 'none')
        xlabel('time t /fs', 'FontSize', fontSize)        
        ylabel('peak current /kA', 'FontSize', fontSize)
%         legend([p1], {['rms mean: ' num2str(mean(x_var), '%5.2f'), ' fs', 10, ...
%                 'FWHM mean: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
%                 'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
%                 'FontSize', fontSize-2)      
              
        set(gca, 'FontSize', fontSize, 'Xlim', lim_x)        

        yyaxis right
        
           p2 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), average_slice_energy_spread_nosase(time_pos_good))
           ylabel('energy spread (Mev)')
        legend([p2(1),p1(1)], { ['energy spread'],['current   ','rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
                        'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
                        'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
                        'FontSize', fontSize-2)
                
        %         

    subplot(2,4,8)
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
        set(gca, 'FontSize', fontSize) 


%%

%% save 

if (flag_save)
        save(['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp2,'_', 'NO_SASE' ,'.mat'], ...
             'comment', 'timestamp', 'name_cam', 'img_bgr_no_sase', 'img_sig_no_sase', 'charge_7FL2XTDS_no_sase', ...
             'timecal_fspixel', 'timecal_fspixel_err', 'timeres_calib', 'ergcal', 'ergcal_err', 'phase_TDS_no_sase', 'tmp2_img', ...
              'x_com', 'x_var', 'x_fwhm', 'x_axis', 'x_profile',...
              'y_com', 'y_var', 'y_fwhm', 'y_axis', 'y_profile');
        %display_message([' - data saved']);
end
%%
    
     % add elog printing button


filename1=['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp1, '.mat']
%text=[comment '\r\n', ...
    %    'data saved in:', filename]


filename2=['/home/ttflinac/user/mflorian/PolariX/FL2/image_calibrated/image_', name_cam, '_', timestamp2,'_', 'NO_SASE' ,'.mat']
text=[comment '\r\n', ...
    'data of LPS with SASE saved in:', filename1 '\r\n', ...
    'data of LPS without SASE saved in:', filename2 ]

     comment='data saved at '
     uicontrol('String', 'Print', 'Callback', @(btn,~)printButtonCallbackDFS(btn, 'FL2 PolariX LPSD@OTR9FL2XTDS', text));

%% Reconstruction



power=zeros( length_x);

power= aver_central_energy_nosase.* aver_current_nosase - aver_central_energy.* aver_current;


figure(1044)
subplot(3,3,2)         
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), power(time_pos_good))

        xlabel('time t /fs', 'FontSize', fontSize)        
        ylabel('peak power (GW)', 'FontSize', fontSize)
end
