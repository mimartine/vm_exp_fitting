clear all
close all
clc

%% Load the data

YCdata = dlmread('YC_Filtered_N4.dat');
ECdata = dlmread('EC_Filtered_N4.dat');
YCdata = [ones(size(YCdata,1),1),YCdata];
ECdata = [2*ones(size(ECdata,1),1),ECdata];
data = [YCdata; ECdata];
data(:,4) = data(:,4) - 1;
data(data(:,3) == 2,:) = [];
data(data(:,3) > 2,3) = data(data(:,3) > 2,3)-1;
asymptote_data = dlmread('asymptote_values.dat');
%
plist = unique(data(:,2));
nsubs = numel(plist);
groups = unique(data(:,4));
ngroups = numel(groups);
PVcol = 6; 

%% Plotting Parameters
[gColor, colorNames] = graphColors(4,0);
tLabels     = {'none-none','rsvp-none','rsvp-rsvp'};
figTitle    = {'Baseline', 'Adaptation', 'Washout', 'Recall'};
parmTitle   = {'Adaptation', 'Washout', 'Recall', 'Savings'};
superLabels = {'YC: None-None','YC:RSVP-None','YC:RSVP-RSVP','EC:None-None','EC:RSVP-None','EC:RSVP-RSVP'};
groupLabels = {'YC','EC'};

lineWidth   = 2;
markerSize  = 1;

xlimit = max(data(:,5));

%%
asymptote_data(:,end+1:end+2) = nan;

for i=1:size(asymptote_data,1)
    baseline_error = data(data(:,2) == asymptote_data(i,2) & data(:,4) == (asymptote_data(i,4) + 1) & data(:,5) >= asymptote_data(i,5),6:7);
    baseline_error = mean(baseline_error,1);
    asymptote_data(i,6:7) = baseline_error;
end % for i... 

mean_baseline = [aggregate(asymptote_data,[1,3,4],6:7,@nanmean),aggregate(asymptote_data,[1,3,4],6:7,@nanstd,1)/sqrt(nsubs/6)];

%%
figure('windowstyle','docked','color','w')

for i = 1:2
    subplot(1,2,i)
    hold on
    errorbar(0.75:1:2.75,...
        mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,3) == 1,4),...
        mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,3) == 1,6), ...
        'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
    errorbar(1.25:1:3.25,...
        mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,3) == 3,4),...
        mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,3) == 3,6), ...        
        'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
    hold off
    
    ylabel('Error (deg)','Fontsize',24)
    xlabel('Condition','Fontsize',24)
    
    set(gca,'XLim', [.5 3+.5], 'XTick', 1:3, 'Fontsize', 18)
    set(gca,'YLim', [0 20],'YTick', 0:5:20);
    set(gca,'XTickLabel',tLabels(1:3),'Fontsize', 12)
    if i == 1
        legend(figTitle([2,4]))
    end
    title(groupLabels{i});

end % for i..

%%
group_mean = [aggregate(data,[1,3,4,5],6:7,@nanmean),aggregate(data,[1,3,4,5],6:7,@nanstd,1)./sqrt(nsubs/6)];

figure('windowStyle','docked','color','w')
cnt = 1;
for i = 1:2
    for j = 1:3
        subplot(2,3,cnt)
        hold on
        errorbar(group_mean(group_mean(:,1) == i & group_mean(:,2) == 1 & group_mean(:,3) == j+1,5),...
            group_mean(group_mean(:,1) == i & group_mean(:,2) == 1 & group_mean(:,3) == j+1,7), ...
            'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
        errorbar(group_mean(group_mean(:,1) == i & group_mean(:,2) == 2 & group_mean(:,3) == j+1,5),...
            group_mean(group_mean(:,1) == i & group_mean(:,2) == 2 & group_mean(:,3) == j+1,7), ...
            'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
        errorbar(group_mean(group_mean(:,1) == i & group_mean(:,2) == 3 & group_mean(:,3) == j+1,5),...
            group_mean(group_mean(:,1) == i & group_mean(:,2) == 3 & group_mean(:,3) == j+1,7), ...
            'ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
        hline = refline([0,mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,2) == 1 & mean_baseline(:,3) == j,4)]);
        hline.Color = gColor(1,:);
        hline.LineWidth = lineWidth;
        hline.LineStyle = '--';
        hline = refline([0,mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,2) == 2 & mean_baseline(:,3) == j,4)]);
        hline.Color = gColor(2,:);
        hline.LineWidth = lineWidth;
        hline.LineStyle = '--';
        hline = refline([0,mean_baseline(mean_baseline(:,1) == i & mean_baseline(:,2) == 3 & mean_baseline(:,3) == j,4)]);
        hline.Color = gColor(3,:);
        hline.LineWidth = lineWidth;
        hline.LineStyle = '--';
        hold off
        
        if j == 1
            set(gca,'XLim', [.5 xlimit + .5], 'XTick', 0:10:xlimit, 'Fontsize', 24)
        else
            set(gca,'XLim', [0.5 xlimit/2 + .5], 'XTick', 0:10:xlimit/2, 'Fontsize',24)
        end
        
        if j == 2
            set(gca,'YLim', [-60 10], 'YTick', -60:10:10, 'Fontsize', 24);
        else
            set(gca,'YLim', [-10 60], 'YTick', -10:10:60, 'Fontsize', 24);
        end
        if cnt == 1
            legend('None-None', 'Low-None' , 'Low-Low');
        end
        cnt = cnt + 1;
    end % for ...
end % for i...

