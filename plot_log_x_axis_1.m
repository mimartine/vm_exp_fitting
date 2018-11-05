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

%
plist = unique(data(:,2));
nsubs = numel(plist);
PVcol = 6; 

%% Plotting Parameters
[gColor, colorNames] = graphColors(4,0);
tLabels     = {'none-none','rsvp-none','rsvp-rsvp'};
figTitle    = {'Baseline', 'Adaptation', 'Washout', 'Recall'};
parmTitle   = {'Adaptation', 'Washout', 'Recall', 'Savings'};
superLabels = {'YC: None-None','YC:RSVP-None','YC:RSVP-RSVP','EC:None-None','EC:RSVP-None','EC:RSVP-RSVP'};

lineWidth   = 2;
markerSize  = 1;
xStart      = [1  41 61];
xEnd        = [40 60 80]; 
xlimit      = [40 20 20];

meanData = [aggregate(data,[1 3 4 5],PVcol,@nanmean),aggregate(data,[1 3 4 5],PVcol,@nanstd,1)./sqrt(nsubs/6)];

%% Plotting
for cIdx = 1:2
    figure('windowstyle','docked','color','w')
    for i = 1:3

        subplot(1,4,i);
        blkData = meanData(meanData(:,1) == cIdx & meanData(:,3) == i+1,:);
        
        x_values = log(1:size(blkData,1)/3);
        
        hold on
        errorbar(x_values,blkData(blkData(:,2) == 1,5), blkData(blkData(:,2) == 1,6), ...
            '-ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
        errorbar(x_values,blkData(blkData(:,2) == 2,5), blkData(blkData(:,2) == 2,6), ...
            '-ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
        errorbar(x_values,blkData(blkData(:,2) == 3,5), blkData(blkData(:,2) == 3,6), ...
            '-ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
        %     errorbar(blkData(blkData(:,1) == 4,4), blkData(blkData(:,1) == 4,5), ...
        %         'ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
        
%         plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 1,(xStart(i):xEnd(i))+2), ...
%             '-ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
%         plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 2,(xStart(i):xEnd(i))+2), ...
%             '-ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
%         plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 3,(xStart(i):xEnd(i))+2), ...
%             '-ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
        %     plot(meanPred(4,(xStart(i):xEnd(i))+1), ...
        %         '-ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
        hold off
        
        ylabel('Error (deg)','Fontsize',24)
        xlabel('Trial Blocks','Fontsize',24)
        
%         set(gca,'XLim', [.5 xlimit(i)+.5], 'XTick', 0:10:xlimit(i), 'Fontsize', 18)
        set(gca,'YLim', [-50 50],'YTick', -50:25:50);
        title(figTitle{i})
        if i == 1
            legend(tLabels);
        end
        
    end
end

%%
for cIdx = 1:2
    for pIdx = 1:3
        figure('windowstyle','docked','color','w')
        cList = unique(data(data(:,1) == cIdx & data(:,3) == pIdx,2));
        for i = 1:length(cList)
            subplot(3,5,i)
            hold on
            d1 = data(data(:,2) == cList(i) & data(:,4) == 2,6);
%             d2 = modelPred(modelPred(:,2) == cList(i),4:43);
            plot(log(1:length(d1)),d1)
%             plot(d2)
            hold off
            set(gca,'YLim', [-10 50],'YTick', -10:10:50);
            title(num2str(cList(i)))
        end % ...
    end
end

%%
% 
% %% Adapatation
% for cIdx = 1:2
%     for pIdx = 1:3
%         figure('windowstyle','docked','color','w')
%         cList = unique(data(data(:,1) == cIdx & data(:,3) == pIdx,2));
%         for i = 1:length(cList)
%             subplot(3,5,i)
%             hold on
%             d1 = data(data(:,2) == cList(i) & data(:,4) == 2,6);
%             d2 = modelPred(modelPred(:,2) == cList(i),4:43);
%             plot(d1)
%             plot(d2)
%             hold off
%             set(gca,'YLim', [-10 50],'YTick', -10:10:50);
%             title(num2str(cList(i)))
%         end % ...
%         supertitle(superLabels{(cIdx-1)*2+pIdx})
%     end
% end
% 
% %% Washout 
% for cIdx = 1:2
%     for pIdx = 1:3
%         figure('windowstyle','docked','color','w')
%         cList = unique(data(data(:,1) == cIdx & data(:,3) == pIdx,2));
%         for i = 1:length(cList)
%             subplot(3,5,i)
%             hold on
%             d1 = data(data(:,2) == cList(i) & data(:,4) == 3,6);
%             d2 = modelPred(modelPred(:,2) == cList(i),44:63);
%             plot(d1)
%             plot(d2)
%             hold off
%             set(gca,'YLim', [-50 10],'YTick', -50:10:10);
%             title(num2str(cList(i)))
%         end % ...
%         supertitle(superLabels{(cIdx-1)*2+pIdx})
%     end
% end
% 
% %% Recall
% for cIdx = 1:2
%     for pIdx = 1:3
%         figure('windowstyle','docked','color','w')
%         cList = unique(data(data(:,1) == cIdx & data(:,3) == pIdx,2));
%         for i = 1:length(cList)
%             subplot(3,5,i)
%             hold on
%             d1 = data(data(:,2) == cList(i) & data(:,4) == 4,6);
%             d2 = modelPred(modelPred(:,2) == cList(i),64:83);
%             plot(d1)
%             plot(d2)
%             hold off
%             set(gca,'YLim', [-10 50],'YTick', -10:10:50);
%             title(num2str(cList(i)))
%         end % ...
%         supertitle(superLabels{(cIdx-1)*2+pIdx})
%     end
% end