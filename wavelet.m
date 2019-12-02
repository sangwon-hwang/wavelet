%% wavelet image extraction

% daily click data %
% dir = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_click_df.csv';
% daily_click = importdata(dir);
% [click_r, click_c] = size(daily_click.data());

dir_5 = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_click_df_5.csv';
daily_click_5 = importdata(dir_5);
[click_r_5, click_c_5] = size(daily_click_5.data());

% figure(1), clf
% hold on
% % subplot
% for i=1:(click_c-1)
%         subplot(2, ceil((click_c-1)/2), i);
%         plot(normalize(daily_click.data(:,i+1)))
%         xlabel('Time (min)'), ylabel('Click')
% end
% hold off

figure(1), clf
hold on
% subplot
for i=1:21
        subplot(3, 7, i);
        plot(normalize(daily_click_5.data(:,i+1)))
        xlabel('Time (min)'), ylabel('Click')
        title([ sprintf('D%d', i) ])
end
hold off


% % raw conversion data %
% dir = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_conversion_df.csv'
% daily_conversion = importdata(dir)
% [conversion_r, conversion_c] = size(daily_conversion.data())

dir_5 = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_conversion_df_5.csv'
daily_conversion_5 = importdata(dir_5)
[conversion_r_5, conversion_c_5] = size(daily_conversion_5.data())

% figure(2), clf
% hold on
% % subplot
% for i=1:(click_c-1)
%     subplot(2, ceil((click_c-1)/2), i);
%     plot(daily_conversion.data(:,i+1))
%     xlabel('Time (min)'), ylabel('Conversion')
% end
% hold off

figure(3), clf
hold on
% subplot
for i=1:21
        subplot(3, 7, i);
        plot(normalize(daily_conversion_5.data(:,i+1)))
        xlabel('Time (min)'), ylabel('Conversion')
        title([ sprintf('D%d', i) ])
end
hold off


% 컨버전 전체
% subplot
%%

for k=1:(conversion_c/5)    
    
    figure(k), clf
    hold on
%     disp(k)
    for i=((k-1)*5+1):(k*5)
%         disp(i)
%         disp((i-(k-1)*5))
        subplot(1, 5, (i-(k-1)*5 ));
        plot(daily_conversion.data(:,i))
        title([ sprintf('%d', i) ])
%         xlabel('Time (min)'), ylabel('Conversion')
    end
    hold off
    saveas(figure(k),sprintf('%d_5day_conversion.jpg', k));

end


% three wavelet transform %
w_type = ["morse", "amor", "bump"]

% magnitude : daily click transform
for k=1:3    
    for j=1:3 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( normalize(daily_click.data(:,j+1), 'range') , w_type(k), 5);
        args = {daily_click.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
        
%         if (j==1|j==11|j==12|j==18)
%             saveas(temp_fig,sprintf('./magnitude/click/train/cat.day%d.click.%s.jpg', j, w_type(k)));
%         elseif (j==20)
%             saveas(temp_fig,sprintf('./magnitude/click/test/cat.day%d.click.%s.jpg', j, w_type(k)));                
%         elseif (13 <= j) && (j <= 19)
%             saveas(temp_fig,sprintf('./magnitude/click/test/dog.day%d.click.%s.jpg', j, w_type(k)));                
%         else
%             saveas(temp_fig,sprintf('./magnitude/click/train/dog.day%d.click.%s.jpg', j, w_type(k)));
%         end        
    end    
end

%% 한개야! temp 5  real(temp_cwt), imag(temp_cwt)
for k=1:3    
    for j=1:1 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_click_5.data(:,j+1) , w_type(3));
        args = {daily_click_5.data(:,1), temp_freq, real(temp_cwt), imag(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,180);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d', j) ])

        
%         if (j==1|j==11|j==12|j==18|j==19|j==20)
%             saveas(temp_fig,sprintf('./magnitude/click/train/cat.day%d.click.%s.jpg', j, w_type(k)));
%         elseif (j==20)
%             saveas(temp_fig,sprintf('./magnitude/click/test/cat.day%d.click.%s.jpg', j, w_type(k)));                
%         elseif (13 <= j) && (j <= 19)
%             saveas(temp_fig,sprintf('./magnitude/click/test/dog.day%d.click.%s.jpg', j, w_type(k)));                
%         else
%             saveas(temp_fig,sprintf('./magnitude/click/train/dog.day%d.click.%s.jpg', j, w_type(k)));
%         end        
    end    
end

% 컨버전

for k=1:3
    for j=1:1 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_conversion_5.data(:,j+1) , w_type(3));
        args = {daily_conversion_5.data(:,1), temp_freq, real(temp_cwt), imag(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,180);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d', j) ])
%         if (j==1|j==11|j==12|j==18)
%             saveas(temp_fig,sprintf('./magnitude/click/train/cat.day%d.click.%s.jpg', j, w_type(k)));
%         elseif (j==20)
%             saveas(temp_fig,sprintf('./magnitude/click/test/cat.day%d.click.%s.jpg', j, w_type(k)));                
%         elseif (13 <= j) && (j <= 19)
%             saveas(temp_fig,sprintf('./magnitude/click/test/dog.day%d.click.%s.jpg', j, w_type(k)));                
%         else
%             saveas(temp_fig,sprintf('./magnitude/click/train/dog.day%d.click.%s.jpg', j, w_type(k)));
%         end        
    end    
end

%%

%% 두개야! temp 5  abs(temp_cwt).^2
% for k=1:3    
    k=2    

    for j=1:21 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_click_5.data(:,j+1) , w_type(k));
        args = {daily_click_5.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d %s', j, w_type(k)) ])

        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/cat.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (2 <= j) && (j <= 10)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==13|j==14)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==19|j==20)
            saveas(temp_fig,sprintf('./magnitude05/%s/test/cat.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        else
            saveas(temp_fig,sprintf('./magnitude05/%s/test/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        end
    end
% end
% temp 5 conversion 5min
% for k=1:3    
    for j=1:21 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_conversion_5.data(:,j+1) , w_type(k));
        args = {daily_conversion_5.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d %s', j, w_type(k)) ])
        
        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/cat.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (2 <= j) && (j <= 10)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==13|j==14)
            saveas(temp_fig,sprintf('./magnitude05/%s/train/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==19|j==20)
            saveas(temp_fig,sprintf('./magnitude05/%s/test/cat.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        else
            saveas(temp_fig,sprintf('./magnitude05/%s/test/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        end
    end    
% end

%% 마지막 세개야! angle(temp_cwt)
% for k=1:3    
k=3
    for j=1:21 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_click_5.data(:,j+1) , w_type(3));
        args = {daily_click_5.data(:,1), temp_freq, angle(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d %s', j, w_type(k)) ])
        
        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./angle05/%s/train/cat.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (2 <= j) && (j <= 10)
            saveas(temp_fig,sprintf('./angle05/%s/train/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==13|j==14)
            saveas(temp_fig,sprintf('./angle05/%s/train/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==19|j==20)
            saveas(temp_fig,sprintf('./angle05/%s/test/cat.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        else
            saveas(temp_fig,sprintf('./angle05/%s/test/dog.day%d.click.%s.jpg', w_type(k), j, w_type(k)));
        end
    end    
% end
% temp 5 conversion 5min
% for k=1:3    
    for j=1:21 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_conversion_5.data(:,j+1) , w_type(3));
        args = {daily_conversion_5.data(:,1), temp_freq, angle(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d %s', j, w_type(k)) ])
        
        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./angle05/%s/train/cat.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (2 <= j) && (j <= 10)
            saveas(temp_fig,sprintf('./angle05/%s/train/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==13|j==14)
            saveas(temp_fig,sprintf('./angle05/%s/train/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        elseif (j==19|j==20)
            saveas(temp_fig,sprintf('./angle05/%s/test/cat.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        else
            saveas(temp_fig,sprintf('./angle05/%s/test/dog.day%d.conv.%s.jpg', w_type(k), j, w_type(k)));
        end
    end    
% end

%%
        figure;
        [temp_cwt, temp_freq] = cwt( daily_click.data(:,2+1) , w_type(3));
        plot(temp_cwt)
        args = {daily_click.data(:,1), temp_freq, real(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';

        figure;
        [temp_cwt, temp_freq] = cwt( daily_conversion.data(:,2+1) , w_type(3));
        args = {daily_click.data(:,1), temp_freq, angle(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';

        
        [wcoh,~,F] = wcoherence(daily_click.data(:,2), daily_conversion.data(:,2), 40);
        figure
        surf(daily_click.data(:,1),F,abs(wcoh).^2)
        shading interp
        axis tight
        hc = colorbar;
        hc.Label.String = 'Coherence';
        title('Wavelet Coherence')
        xlabel('Seconds')
        ylabel('Hz')
        ylim([0 2.5])
        set(gca,'ytick',[0.15 1.2 2])




%  real : daily conversion transform
for k=1:3    
    for j=1:20
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt(daily_conversion.data(:,j+1), w_type(k));
        args = {daily_click.data(:,1), temp_freq, real(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
        if (1 <= j) && (j <= 13)
        saveas(temp_fig,sprintf('./real/train/dog.day%d.conv.%s.jpg', j, w_type(k)));
        else
        saveas(temp_fig,sprintf('./real/test/dog.day%d.conv.%s.jpg', j, w_type(k)));            
        end
    end    
end


    
    
    
    
% subplot all together : magnitude : daily click transform
for k=1:3    
    temp_fig = figure;
    for j=1:(click_c-1) % 7 days 
        [temp_cwt, temp_freq] = cwt( normalize(daily_click.data(:,j+1), 'range') , w_type(k));
        subplot(1, ceil((click_c-1)/2), j);
        args = {daily_click.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
    end
    saveas(temp_fig,'mag.daily.click.%s.jpg', j, w_type(k));
end

% subplot all together : magnitude : daily conversion transform
for k=1:3    
    temp_fig = figure;
    for j=1:(conversion_c-1)  
        [temp_cwt, temp_freq] = cwt( normalize(daily_conversion.data(:,j+1), 'range') , w_type(k));
        subplot(1, ceil((click_c-1)/2), j);
        args = {daily_click.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
    end
    saveas(temp_fig,'mag.daily.conversion.%s.jpg', j, w_type(k));
end

% subplot all together : real : daily click transform
for k=1:3    
    temp_fig = figure;
    for j=1:(click_c-1) % 7 days 
        [temp_cwt, temp_freq] = cwt( normalize(daily_click.data(:,j+1), 'range') , w_type(k));
        subplot(1, ceil((click_c-1)/2), j);
        args = {daily_click.data(:,1), temp_freq, real(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
    end
    saveas(temp_fig,'real.daily.click.%s.jpg', j, w_type(k));
end

% subplot all together : magnitude : daily conversion transform
for k=1:3    
    temp_fig = figure;
    for j=1:(conversion_c-1)  
        [temp_cwt, temp_freq] = cwt( normalize(daily_conversion.data(:,j+1), 'range') , w_type(k));
        subplot(1, ceil((click_c-1)/2), j);
        args = {daily_click.data(:,1), temp_freq, real(temp_cwt)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
    end
    saveas(temp_fig,'real.daily.conversion.%s.jpg', j, w_type(k));
end


        
%%






figure
[temp_cwt, temp_freq] = cwt(normalize(daily_click.data(:,5)), 2e4, 'VoicesPerOctave', 48);
        subplot(1,1,1);
        args = {daily_click.data(:,1), temp_freq, abs(temp_cwt).^2};
%         surf(args{:},'edgecolor','none');
        contourf(args{:},'edgecolor','none');
%         view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'click';


figure
[temp_cwt, temp_freq] = cwt(normalize(daily_conversion.data(:,5)), 2e4, 'VoicesPerOctave',  48);
        subplot(1,1,1);
        args = {daily_conversion.data(:,1), temp_freq, abs(temp_cwt).^2};
%         surf(args{:},'edgecolor','none');
        contourf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'conversion';


        
        
[cwt_click(i),click_f(i)] = cwt(daily_click.data(:,i+1),'morse');
% last number is frequency
% second parameter is type of wavelet
% 'morse', 'amor', and 'bump' - Morse(default), Morlet (Gabor), and bump wavelet

% subplot(7, 3, );
args = {daily_click.data(:,1), click_01_f, abs(cwt_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';

[CWT_click_01,click_01_f] = cwt(daily_click.data(:,2),'amor');
subplot(1, 3, 2);
args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';

[CWT_click_01,click_01_f] = cwt(daily_click.data(:,2),'bump');
subplot(1, 3, 3);
args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';
% end 

% day2 
figure(7), clf
[CWT_click_01,click_01_f] = cwt(daily_click.data(:,3),'morse');
% last number is frequency
% second parameter is type of wavelet
% 'morse', 'amor', and 'bump' - Morse(default), Morlet (Gabor), and bump wavelet
subplot(1, 3, 1);
args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';

[CWT_click_01,click_01_f] = cwt(daily_click.data(:,3),'amor');
subplot(1, 3, 2);
args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';

[CWT_click_01,click_01_f] = cwt(daily_click.data(:,3),'bump');
subplot(1, 3, 3);
args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';













figure(11), clf
[CWT_click_01,click_01_f] = cwt(daily_click.data(:,2),'morse',2e4,'VoicesPerOctave',16);

args = {daily_click.data(:,1),click_01_f,abs(CWT_click_01).^2};
% options : imag(CWT_click_01), real(CWT_click_01) 
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
shading interp; colormap(parula(128));
h = colorbar;
h.Label.String = 'Power';

figure(4), clf
[dpoaeCWT,f] = cwt(dpoaets,2e4,'VoicesPerOctave',16);
helperCWTTimeFreqPlot(dpoaeCWT,t.*1000,f,...
    'surf','CWT of OAE','milliseconds','Hz')


%%
% function helperCWTTimeFreqPlot(cfs,time,freq,PlotType,varargin)
%   This function helperCWTTimeFreqPlot is only in support of
%   CWTTimeFrequencyExample and PhysiologicSignalAnalysisExample. 
%   It may change in a future release.

params = parseinputs(varargin{:});


%     if strncmpi(PlotType,'surf',1)
%         args = {time,freq,abs(cfs).^2};
%         surf(args{:},'edgecolor','none');
%         view(0,90);
%         axis tight;
%         shading interp; colormap(parula(128));
%         h = colorbar;
%         h.Label.String = 'Power';
%             if isempty(params.xlab) && isempty(params.ylab)
%                 xlabel('Time'); ylabel('Hz');
%             else
%              xlabel(params.xlab); ylabel(params.ylab);
%             end
% 



figure
[CWT_click_day_two,f] = cwt(daily_click.data(:,2),2e4,'VoicesPerOctave',16);
% helperCWTTimeFreqPlot(CWT_click_day_one,f,'surf','EEG','Time','Hz')
contourf(real(CWT_click_day_two),'linecolor','k')

figure
[CWT_click_day_three,f] = cwt(daily_click.data(:,3),2e4,'VoicesPerOctave',16);
% helperCWTTimeFreqPlot(CWT_click_day_one,f,'surf','EEG','Time','Hz')
contourf(real(CWT_click_day_three),'linecolor','k')

figure
[CWT_click_day_four,f] = cwt(daily_click.data(:,4),2e4,'VoicesPerOctave',16);
% helperCWTTimeFreqPlot(CWT_click_day_one,f,'surf','EEG','Time','Hz')
contourf(real(CWT_click_day_four),'linecolor','k')



figure
plot(daily_conversion.data(:,2))

figure
[CWT_conversion_day_one,f] = cwt(daily_conversion.data(:,2),2e4,'VoicesPerOctave',16);
% helperCWTTimeFreqPlot(CWT_click_day_one,f,'surf','EEG','Time','Hz')
contourf(real(CWT_conversion_day_one),'linecolor','k')


figure
[CWT_conversion_day_two,f] = cwt(daily_conversion.data(:,3),2e4,'VoicesPerOctave',16);
% helperCWTTimeFreqPlot(CWT_click_day_one,f,'surf','EEG','Time','Hz')
contourf(real(CWT_conversion_day_two),'linecolor','k')


figure
[dataCWT,f] = cwt(data,2e4,'VoicesPerOctave',16);
helperCWTTimeFreqPlot(dataCWT,f,'surf','EEG','Time','Hz')



dailyconversiondf.data






%%
% You can investigate the time evolution of the OAE by finding the CWT
% coefficients closest in frequency to 1230 Hz and examining their
% magnitudes as a function of time. Plot the magnitudes along with time
% markers designating the beginning and end of the evoking stimulus.

[~,idx1230] = min(abs(f-1230));
cfsOAE = dpoaeCWT(idx1230,:);
plot(t.*1000,abs(cfsOAE))
hold on
AX = gca;
plot([25 25],[AX.YLim(1) AX.YLim(2)],'r')
plot([175 175],[AX.YLim(1) AX.YLim(2)],'r')
xlabel('msec')
title('CWT Coefficient Magnitudes')
%%
% There is some delay between the onset of the evoking stimulus and the
% OAE. Once the evoking stimulus is terminated, the OAE immediately begins
% to decay in magnitude.
%%
% Another way to isolate the emission is to use the inverse CWT to
% reconstruct a frequency-localized approximation in the time domain.
%
% Construct a frequency-localized emission approximation by extracting the
% CWT coefficients corresponding to frequencies between 1150 and 1350 Hz.
% Use these coefficients and invert the CWT. Plot the original data along
% with the reconstruction and markers indicating the beginning and end of
% the evoking stimulus.

frange = [1150 1350];
xrec = icwt(dpoaeCWT,f,frange);
figure
plot(t.*1000,dpoaets)
hold on
plot(t.*1000,xrec,'r')
AX = gca;
ylimits = AX.YLim;
plot([25 25],ylimits,'k')
plot([175 175],ylimits,'k')
grid on
xlabel('Milliseconds')
ylabel('Amplitude')
title('Frequency-Localized Reconstruction of Emission')
%%
% In the time-domain data, you clearly see how the emission ramps on and
% off at the application and termination of the evoking stimulus.
%
% It is important to note that even though a range of frequencies were
% selected for the reconstruction, the analytic wavelet transform actually
% encodes the exact frequency of the emission. To demonstrate this, take
% the Fourier transform of the emission approximation reconstructed from
% the analytic CWT.

xdft = fft(xrec);
freq = 0:2e4/numel(xrec):1e4;
xdft = xdft(1:numel(xrec)/2+1);
figure
plot(freq,abs(xdft))
xlabel('Hz')
ylabel('Magnitude')
title('Fourier Transform of CWT-Based Signal Approximation')
[~,maxidx] = max(abs(xdft));
fprintf('The frequency is %4.2f Hz\n',freq(maxidx))
%%
% This example used |cwt| to obtain a time-frequency analysis of the OAE
% data and <matlab:doc('icwt') icwt> to obtain a frequency-localized
% approximation to the signal.


%% References
% 
% Cui, X., Bryant, D.M., and Reiss. A.L. "NIRS-Based hyperscanning reveals
% increased interpersonal coherence in superior frontal cortex during
% cooperation", Neuroimage, 59(3), 2430-2437, 2012.
%
% Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
% Mietus JE, Moody GB, Peng C-K, Stanley HE. "PhysioBank, PhysioToolkit,
% and PhysioNet: Components of a New Research Resource for Complex
% Physiologic Signals." Circulation 101(23):e215-e220, 2000. 
% |http://circ.ahajournals.org/cgi/content/full/101/23/e215|
%   
% Mallat, S. "A Wavelet Tour of Signal Processing: The Sparse Way",
% Academic Press, 2009.
%
% Moody, G.B. "Evaluating ECG Analyzers".
% |http://www.physionet.org/physiotools/wfdb/doc/wag-src/eval0.tex|
%
% Moody GB, Mark RG. "The impact of the MIT-BIH Arrhythmia Database." IEEE
% Eng in Med and Biol 20(3):45-50 (May-June 2001).
%% Appendix
% The following helper functions are used in this example.
%
% * <matlab:edit('helperCWTTimeFreqPlot.m') helperCWTTimeFreqPlot.m>
% * <matlab:edit('helperCWTTimeFreqVector.m') helperCWTTimeFreqVector.m>
