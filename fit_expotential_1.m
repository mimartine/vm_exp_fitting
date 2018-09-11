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

%% Fitting Parameters
parms = -0.25;

minParms                    = -1;
maxParms                    = 1;
parmSteps                   = 0.1;
psrFactor                   = 0.5;
psaFactor                   = 1.5;

%
nTrials     = sum(xlimit);
modelPred   = nan(nsubs,nTrials);

model = @fit_expcurve;

%% Fitting
for pIdx = 1:nsubs
    for cIdx = 1:3
        pData = data(data(:,2) == plist(pIdx) & data(:,4) == cIdx+1,:);
        pno = mean(pData(:,2));
        cnd = mean(pData(:,3));
        grp = mean(pData(:,1));
              
        pData   = pData(:,PVcol)';
        nt      = length(pData);
        trials  = 1:nt; 
        trials  = trials(~isnan(pData));
        
        pData(isnan(pData)) = [];
        
%         if cIdx ~= 2
%             peakPt = max(pData);
%         else
%             peakPt = min(pData);
%         end
%         
        peakPt = pData(1);
      
        asymPt  = mean(pData(end-3:end));
        delta   = peakPt - asymPt;
        
        %     options = optimset('Display','final','MaxIter',10000,'MaxFunEvals',1e16,'TolFun',1e-16,'TolX',1e-16);
        %     [fit, fval] = fmincon(model,parms,Aeq,beq,[],[],minParms,maxParms,[],options, [sp delta], pData);
        
        [fit, fval, exitflag] =...
            hkjeeves(model, parms,...
            0, minParms, maxParms,...  % Set first parm to -1 to show iteration
            -1, 1e-4, 10000,...
            parmSteps, psaFactor, psrFactor,...
            [asymPt delta], pData,trials);
        
        fValue(pIdx,cIdx) = fval;
        fParms(pIdx,cIdx) = fit;
        
        y = expFun([asymPt delta fit],nt);
        
        r = power(corrcoef([pData;y(trials)]'),2);
        rSquared(pIdx,cIdx) = r(1,2);
%         
%         if length(pData) < nt
%             y = expFun([asymPt delta fit],nt);
%         end
              
        modelPred(pIdx,xStart(cIdx):xEnd(cIdx)) = y;
    end % for j...
end % for i...

%%
modelPred = [aggregate(data,1:3,1:3,@mean,1), modelPred];

meanPred = aggregate(modelPred,[1 3],4:size(modelPred,2),@nanmean);
meanData = [aggregate(data,[1 3 4 5],PVcol,@nanmean),aggregate(data,[1 3 4 5],PVcol,@nanstd,1)./sqrt(nsubs/6)];

fParms = [aggregate(data,1:3, 1:3,@mean,1),-(fParms)];
fParms(:,7) = fParms(:,6) - fParms(:,4);
meanParms = aggregate(fParms,[1 3],4:7);
stdParms  = aggregate(fParms,[1 3],4:7,@std);
stdParms(:,4:end) = stdParms(:,4:end)./sqrt(nsubs/6);

%% Plotting


for cIdx = 1:2
    figure('windowstyle','docked','color','w')
    for i = 1:3
        subplot(1,4,i);
        blkData = meanData(meanData(:,1) == cIdx & meanData(:,3) == i+1,:);
        
        hold on
        errorbar(blkData(blkData(:,2) == 1,5), blkData(blkData(:,2) == 1,6), ...
            'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
        errorbar(blkData(blkData(:,2) == 2,5), blkData(blkData(:,2) == 2,6), ...
            'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
        errorbar(blkData(blkData(:,2) == 3,5), blkData(blkData(:,2) == 3,6), ...
            'ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
        %     errorbar(blkData(blkData(:,1) == 4,4), blkData(blkData(:,1) == 4,5), ...
        %         'ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
        
        plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 1,(xStart(i):xEnd(i))+2), ...
            '-ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
        plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 2,(xStart(i):xEnd(i))+2), ...
            '-ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
        plot(meanPred(meanPred(:,1) == cIdx & meanPred(:,2) == 3,(xStart(i):xEnd(i))+2), ...
            '-ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
        %     plot(meanPred(4,(xStart(i):xEnd(i))+1), ...
        %         '-ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
        hold off
        
        ylabel('Error (deg)','Fontsize',24)
        xlabel('Trial Blocks','Fontsize',24)
        
        set(gca,'XLim', [.5 xlimit(i)+.5], 'XTick', 0:10:xlimit(i), 'Fontsize', 18)
        set(gca,'YLim', [-50 50],'YTick', -50:25:50);
        title(figTitle{i})
        if i == 1
            legend(tLabels);
        end
        
    end
end

%%
for cIdx = 1:2
    figure('windowstyle','docked','color','w')
    for i = 1:4
        subplot(2,2,i);
        hold on
        errorbar(1,meanParms(meanParms(:,1) == cIdx & meanParms(:,2) == 1,i+2),...
            stdParms(stdParms(:,1) == cIdx & stdParms(:,2) == 1,i+2),...
            'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize*4,'MarkerFaceColor',gColor(1,:));
        errorbar(2,meanParms(meanParms(:,1) == cIdx & meanParms(:,2) == 2,i+2),...
            stdParms(stdParms(:,1) == cIdx & stdParms(:,2) == 2,i+2),...
            'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize*4,'MarkerFaceColor',gColor(2,:));
        errorbar(3,meanParms(meanParms(:,1) == cIdx & meanParms(:,2) == 3,i+2),...
            stdParms(stdParms(:,1) == cIdx & stdParms(:,2) == 3,i+2),...
            'ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize*4,'MarkerFaceColor',gColor(3,:));
        hold off
        
        ylabel('Learning Rate','Fontsize',24)
        xlabel('Group','Fontsize',24)
        
        set(gca,'XLim', [.5 3.5], 'XTick', 1:3, 'XTickLabel', tLabels, 'Fontsize', 16)
        
        if i ~= 4
            set(gca,'YLim', [0 1],'YTick', 0.0:0.25:1,'Fontsize', 16);
        else
            set(gca,'YLim', [-0.5 0.5],'YTick', -0.5:0.25:0.5,'Fontsize', 16);
        end
        title(parmTitle{i})
    end % for i...
end

for cIdx = 1:2
    for pIdx = 1:3
        figure('windowstyle','docked','color','w')
        cList = unique(data(data(:,1) == cIdx & data(:,3) == pIdx,2));
        for i = 1:length(cList)
            subplot(3,5,i)
            hold on
            d1 = data(data(:,2) == cList(i) & data(:,4) == 2,6);
            d2 = modelPred(modelPred(:,2) == cList(i),4:43);
            plot(d1)
            plot(d2)
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