
% pulse reconstrcution function 
% Najmeh mirian 29 Nov 2023

function filename= Pulse_reconstruction(a)

% Bakg=load('background_img.mat');
imag_lasing_on=load('imag_shots_lasing_on.mat');
imag_lasing_off=load('imag_shots_lasing_off.mat');

% calibration data
load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'timecal_fspixel', ...
    'timecal_fspixel_err', 'timecal_streak', 'time_res');
load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err')
display('data is loaded');

%num_sig=app.num_shot

%img_sig    = imag_file.img(end:a, :, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% analysis
%img_bgr_mean=Bakg.img_bgr_mean;
% img_sig=imag_lasing_on.img;
% charge_7FL2XTDS_mean=mean(mag_lasing_on.charge.data)  ; %V check it
fontSize        = 14;
calib_x             = timecal_fspixel;


calib_y             = ergcal;

timeres_calib       = time_res;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for on_off=0:3:3
    if on_off==0
        img_sig=imag_lasing_on.img(end:a, :, :);
        num_sig     = size(img_sig, 3);
        for i=1:num_sig, charge_7FL2XTDS(i)=imag_lasing_on.charge(i).data; end
        charge_7FL2XTDS_mean=mean(charge_7FL2XTDS)  ; %V check it

        text='lasing on'
    else
        img_sig=imag_lasing_off.img(end:a, :, :);
        %         charge_7FL2XTDS_mean=mean(imag_lasing_off.charge)  ; %V check it
        %         charge_7FL2XTDS=imag_lasing_off.charge;
        for i=1:num_sig, charge_7FL2XTDS(i)=imag_lasing_off.charge(i).data; end
        charge_7FL2XTDS_mean=mean(charge_7FL2XTDS)  ; %V check it
        text='lasing off'

    end


    img_filt    = img_sig;

    num_sig     = size(img_sig, 3);
    length_y    = size(img_sig, 1);
    length_x    = size(img_sig, 2);

    %img_bgr_mean = squeeze(mean(img_bgr, 3));

    x_com       = zeros([1, num_sig]);              y_com       = zeros([1, num_sig]);
    x_var       = zeros([1, num_sig]);              y_var       = zeros([1, num_sig]);
    x_fwhm      = zeros([1, num_sig]);              y_fwhm      = zeros([1, num_sig]);
    x_profile   = zeros([num_sig, length_x]);       y_profile   = zeros([num_sig, length_y]);
    y=1:length_y;
    %%
    for jj = 1:num_sig

        tmp_img             = hlc_clean_image( squeeze(img_sig(:,:,jj)) );
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
    tmp2_img=squeeze(transpose(img_filt(:,:,end)));

    %% plot 1

    erg_axis            = ergcal * ((1:length_y)-length_y/2);
    time_axis           = timecal_fspixel * ((1:length_x)-length_x/2);
    erg_pos_good =  find(y_profile(end,:)>0); erg_pos_good = min(erg_pos_good):max(erg_pos_good);
    time_pos_good = find(x_profile(end,:)>0); time_pos_good = min(time_pos_good):max(time_pos_good);
    cent_energy_cr=zeros(num_sig, length_x);
    slice_enrgy_spread=zeros(num_sig, length_x);
    %%
    %     figure; plot(slice_enrgy_spread_ave)
    %%
    for jj=1:num_sig
        for n=min(time_pos_good):max(time_pos_good)

          %  mu0=  find( img_filt(erg_pos_good,n,jj)==max(img_filt(erg_pos_good,n,jj)));
          %  [sigma, mu] = gaussfit( erg_pos_good, smoothdata(img_filt(erg_pos_good,n,jj)) , 0, erg_pos_good(mu0(1)));
          [baseline, height,sigma, mu]=hlc_fit_gaussian(x, y);
          if isnan(sigma)
            else

                slice_enrgy_spread(jj, n)=sigma*calib_y;
                cent_energy_cr(jj, n)=mu*calib_y;
            end

        end
    end


    %% jetting of stabilitty correction
    for jj=1: num_sig
        ftemp=find(x_profile(jj,:)==max(x_profile(jj,:)));
        mass_center(jj)=ftemp(1);

    end
    %%
    figure(1044);hold on
    set(gcf, 'position', [100,1, 1500, 1000])
    sumes=zeros( 1,length_x);
    sumEcr=zeros( 1,length_x);
    sumCu=zeros( 1,length_x);

    for jj=1:num_sig
        subplot(3,3,1+on_off);hold on
        plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),slice_enrgy_spread(jj,time_pos_good-mass_center(1)+mass_center(jj)))
        sumes(time_pos_good)=sumes(time_pos_good)+slice_enrgy_spread(jj,time_pos_good-mass_center(1)+mass_center(jj));
        subplot(3,3,2+on_off);hold on
        plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),cent_energy_cr(jj,time_pos_good-mass_center(1)+mass_center(jj)))
        sumEcr(time_pos_good)=sumEcr(time_pos_good)+cent_energy_cr(jj,time_pos_good-mass_center(1)+mass_center(jj));
        subplot(3,3,3+on_off);hold on
        plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),1e-3*charge_7FL2XTDS(jj)*1e-9*1e15*x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj)))
        sumCu(time_pos_good)=sumCu(time_pos_good)+1e-3*charge_7FL2XTDS(jj)*1e-9*1e15*x_profile(jj,time_pos_good-mass_center(1)+mass_center(jj));
    end
    aver_slice_energy_spread=sumes/num_sig;
    aver_central_energy=sumEcr/num_sig;
    aver_current=sumCu/num_sig;

    subplot(3,3,1+on_off)
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_slice_energy_spread(time_pos_good),'--k', 'LineWidth',3 )

    set(gca, 'FontSize', fontSize,'YDir', 'normal')
    xlabel('time t /fs', 'FontSize', fontSize)
    ylabel('slice energy Sepread (Mev)', 'FontSize', fontSize)


    subplot(3,3,2+on_off)
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_central_energy(time_pos_good),'--k', 'LineWidth',3 )

    set(gca, 'FontSize', fontSize,'YDir', 'normal')
    xlabel('time t /fs', 'FontSize', fontSize)
    ylabel('\Delta energy cneter (Mev)', 'FontSize', fontSize)
    title(text, 'FontSize', fontSize+3)
    subplot(3,3,3+on_off)
    plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_current(time_pos_good),'--k', 'LineWidth',3 )

    set(gca, 'FontSize', fontSize,'YDir', 'normal')
    lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');
    xlabel('time t /fs', 'FontSize', fontSize)
    ylabel('peak current /kA', 'FontSize', fontSize)
    %%
    if on_off==0
        imag_lasing_on.aver.central_energy=aver_central_energy;
        imag_lasing_on.aver.current =aver_current;
        imag_lasing_on.aver.slice_energy_spread=aver_slice_energy_spread;
    else
        on_off=4;
        imag_lasing_off.aver.central_energy=aver_central_energy;
        imag_lasing_off.aver.current =aver_current;
        imag_lasing_off.aver.slice_energy_spread=aver_slice_energy_spread;

    end

    % disp(' - ')
    figure(1)
    set(gcf, 'Position', [800, 1, 1500, 1600])

    subplot(2,4,1+on_off)
    imagesc(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'),...
        erg_axis(erg_pos_good), flipud(tmp2_img(time_pos_good,erg_pos_good)') )
    grid on
    title(['single image ', text], 'FontSize', fontSize, 'Interpreter', 'none')
    xlabel('time t /fs)', 'FontSize', fontSize)
    ylabel('rel. energy deviation \delta / 10^{-3}' , 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize,'YDir', 'normal')
    lim_y = get(gca, 'YLim'); lim_x = get(gca, 'XLim');

    subplot(2,4,2+on_off)
    p1 = plot(y_profile(end,erg_pos_good), -erg_axis(erg_pos_good) );
    grid on
    title(text, 'FontSize', fontSize)
    xlabel('norm. erg density', 'FontSize', fontSize)
    ylabel('rel. energy deviation \delta / 10^{-3}', 'FontSize', fontSize)
    legend([p1(1)], {['rms: ' num2str(y_var(end), '%5.2f'), ' \cdot 10^{-3}', 10, ...
        'FWHM: ', num2str(y_fwhm(end), '%5.2f'), ' \cdot 10^{-3}']}, ...
        'Location','northoutside','FontSize', fontSize-2)
    set(gca, 'FontSize', fontSize, 'YDir', 'normal');%, 'Ylim', -lim_y)

    subplot(2,4,3+on_off)

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
        'Location','northoutside','FontSize', fontSize-2)
    yyaxis right

    p2 = plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), aver_slice_energy_spread(time_pos_good))

    ylabel('energy spread (Mev)')
    legend([p2(1),p1(1)], { ['energy spread'],['current   ','rms: ' num2str(x_var(end), '%5.2f'), ' fs', 10, ...
        'FWHM: ', num2str(x_fwhm(end), '%5.2f'), ' fs', 10, ...
        'Time res: ' num2str((timeres_calib), '%5.1f'), ' fs']}, ...
        'FontSize', fontSize-2, 'Location','northoutside')


    set(gca, 'FontSize', fontSize, 'Xlim', lim_x)
    %'FWHM: ', num2str(mean(x_fwhm), '%5.2f'), ' fs', 10, ...
    %

    subplot(2,4,4+on_off)
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
        'Location','northoutside', 'FontSize', fontSize-2)

    set(gca, 'FontSize', fontSize)




end

%% Reconstruction


mass_center_on=find(imag_lasing_on.aver.current==max(imag_lasing_on.aver.current));
mass_center_off=find(imag_lasing_off.aver.current==max(imag_lasing_off.aver.current));


power=zeros( length_x,1);

power= imag_lasing_off.aver.central_energy(time_pos_good).* imag_lasing_off.aver.current(time_pos_good) - ...
    imag_lasing_on.aver.central_energy(time_pos_good-mass_center_off+mass_center_on).* imag_lasing_on.aver.current(time_pos_good-mass_center_off+mass_center_on);
% power=aver_central_energy_lasingoff(time_pos_good).* aver_current_lasingoff(time_pos_good) - ...
%         aver_central_energy_lasingon(time_pos_good).* aver_current_lasingoff(time_pos_good);
power2=(imag_lasing_off.aver.central_energy(time_pos_good)- ...
    imag_lasing_on.aver.central_energy(time_pos_good-mass_center_off+mass_center_on)).* imag_lasing_on.aver.current(time_pos_good-mass_center_off+mass_center_on);
%% save
timeCallib=load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat');
energyCalib=load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat');
t=date;
Year=year(t);
fulder_add=['/home/ttflinac/data/PolariX/', num2str(Year), '/'];
timestamp       = datestr(clock, 'yyyy-mm-ddTHHMMSS');
save([fulder_add, timestamp, '.mat'],'timestamp' ,  'img_bgr_mean',  'imag_lasing_on', 'imag_lasing_off','energyCalib', 'timeCallib', 'power')

display(' - data saved');
%%%%%%%%%%%%%%%%%%%%%
%%


figure(1044)
subplot(3,3,7)
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), power(:), '-k', 'LineWidth',3 )
hold on
plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), power2(:), '-b', 'LineWidth',3 )
xlabel('time t /fs', 'FontSize', fontSize)
ylabel('peak power (GW)', 'FontSize', fontSize)
set(gca, 'FontSize', fontSize,'YDir', 'normal')

subplot(3,3,8)
hold on
p1=plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), imag_lasing_off.aver.central_energy(time_pos_good), '--')
p2=plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), imag_lasing_on.aver.central_energy(time_pos_good-mass_center_off+mass_center_on), '--')

set(gca, 'FontSize', fontSize,'YDir', 'normal')
xlabel('time t /fs', 'FontSize', fontSize)
ylabel('\Delta energy cneter (Mev)', 'FontSize', fontSize)

legend([p1(1), p2(1)],{'lassing off', 'lasing on'}, 'FontSize', fontSize-2, 'Location','northoutside', 'Orientation','horizontal')
figure(1044)
subplot(3,3,9)

hold on
p1=plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), imag_lasing_off.aver.current(time_pos_good), '--')
p2=plot(time_axis(time_pos_good)-mean(time_axis(time_pos_good),'omitnan'), imag_lasing_on.aver.current(time_pos_good-mass_center_off+mass_center_on), '--')
xlabel('time t /fs', 'FontSize', fontSize)
ylabel('Current(KA)', 'FontSize', fontSize)
set(gca, 'FontSize', fontSize,'YDir', 'normal')
legend([p1(1), p2(1)],{'lassing off', 'lasing on'},'FontSize', fontSize-2, 'Location','northoutside', 'Orientation','horizontal')


filename=[fulder_add, timestamp, '.mat']