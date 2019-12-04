%% wavelet image extraction
% daily click data where the time unit is 1 min
dir = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_click_df.csv';
daily_click = importdata(dir);
[click_r, click_c] = size(daily_click.data());

% daily click data where the time unit is 5 min
dir_5 = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_click_df_5.csv';
daily_click_5 = importdata(dir_5);
[click_r_5, click_c_5] = size(daily_click_5.data());

figure(1), clf
hold on
% subplot
for i=1:(click_c-1)
        subplot(2, ceil((click_c-1)/2), i);
        plot(normalize(daily_click.data(:,i+1)))
        xlabel('Time (min)'), ylabel('Click')
end
hold off

figure(2), clf
hold on
% subplot
for i=1:21
        subplot(3, 7, i);
        plot(normalize(daily_click_5.data(:,i+1)))
        xlabel('Time (min)'), ylabel('Click')
        title([ sprintf('D%d', i) ])
end
hold off

%%
% daily conversion data where the time unit is 1 min
dir = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_conversion_df.csv'
daily_conversion = importdata(dir)
[conversion_r, conversion_c] = size(daily_conversion.data())

% daily conversion data where the time unit is 5 min
dir_5 = '/Users/sangwonhwang/Desktop/mikecohen_퓨리에변환/data/daily_conversion_df_5.csv'
daily_conversion_5 = importdata(dir_5)
[conversion_r_5, conversion_c_5] = size(daily_conversion_5.data())

figure(3), clf % each figure has 5 sub plots
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

figure(4), clf % each figure has 5 sub plots
for k=1:(conversion_c5/5)    
    
    figure(k), clf
    hold on
%     disp(k)
    for i=((k-1)*5+1):(k*5)
%         disp(i)
%         disp((i-(k-1)*5))
        subplot(1, 5, (i-(k-1)*5 ));
        plot(daily_conversion_5.data(:,i))
        title([ sprintf('%d', i) ])
%         xlabel('Time (min)'), ylabel('Conversion')
    end
    hold off
    saveas(figure(k),sprintf('%d_5day_conversion.jpg', k));

end

% figure(3), clf
% hold on
% % subplot
% for i=1:(click_c-1)
%     subplot(2, ceil((click_c-1)/2), i);
%     plot(daily_conversion.data(:,i+1))
%     xlabel('Time (min)'), ylabel('Conversion')
% end
% hold off
% 
% figure(4), clf
% hold on
% % subplot
% for i=1:21
%         subplot(3, 7, i);
%         plot(normalize(daily_conversion_5.data(:,i+1)))
%         xlabel('Time (min)'), ylabel('Conversion')
%         title([ sprintf('D%d', i) ])
% end
% hold off

%%
% three wavelet transform %
w_type = ["morse", "amor", "bump"]

%%
% magnitude : daily click transform where the time unit is 1 min
for k=1:3    
    for j=1:21 % 3 weeks 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( normalize(daily_click.data(:,j+1), 'range') , w_type(k), 5);
        args = {daily_click.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(parula(128));
        h = colorbar;
        h.Label.String = 'Power';
        
        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./magnitude/click/train/cat.day%d.click.%s.jpg', j, w_type(k)));
        elseif (j==20)
            saveas(temp_fig,sprintf('./magnitude/click/test/cat.day%d.click.%s.jpg', j, w_type(k)));                
        elseif (13 <= j) && (j <= 19)
            saveas(temp_fig,sprintf('./magnitude/click/test/dog.day%d.click.%s.jpg', j, w_type(k)));                
        else
            saveas(temp_fig,sprintf('./magnitude/click/train/dog.day%d.click.%s.jpg', j, w_type(k)));
        end        
    end    
end

%%
% magnitude : daily click transform where the time unit is 5 mins
for k=1:3    
    for j=1:7 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_click_5.data(:,j+1) , w_type(k));
        args = {daily_click_5.data(:,1), temp_freq, abs(temp_cwt).^2};
        surf(args{:},'edgecolor','none');
        view(0,180);
        axis tight;
        shading interp; colormap(jet(128));
        h = colorbar;
        h.Label.String = 'Power';
        title([ sprintf('D%d', j) ])
        
        if (j==1|j==11|j==12|j==18)
            saveas(temp_fig,sprintf('./magnitude/click/train/cat.day%d.click.%s.jpg', j, w_type(k)));
        elseif (j==20)
            saveas(temp_fig,sprintf('./magnitude/click/test/cat.day%d.click.%s.jpg', j, w_type(k)));                
        elseif (13 <= j) && (j <= 19)
            saveas(temp_fig,sprintf('./magnitude/click/test/dog.day%d.click.%s.jpg', j, w_type(k)));                
        else
            saveas(temp_fig,sprintf('./magnitude/click/train/dog.day%d.click.%s.jpg', j, w_type(k)));
        end        
    end
end    

%% 
% magnitude : daily click transform where the time unit is 5 mins with different if-else
for k=1:3    
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
end

%% 
% angle : daily click transform where the time unit is 5 mins
for k=1:3    
    for j=1:21 % 7 days 
        temp_fig = figure;
        [temp_cwt, temp_freq] = cwt( daily_click_5.data(:,j+1) , w_type(k));
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
end

%%
% 3 dimensional daily conversion transform : time, x-axis(real value), y-axis(imaginary value)
%                                            where the time unit is 5 mins
for k=1:3    
    for j=1:7 % 7 days 
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
% magnitude : daily conversion transform where the time unit is 5 mins
for k=1:3    
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
end

%%
% angle : daily conversion transform where the time unit is 5 mins
for k=1:3    
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
end

%% wavelet coefficience test      

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
%%