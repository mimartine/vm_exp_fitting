%% Fit Expotential Function for Experiment 1 of Prediction Error experiment
% Single parameter model -> rate of adaptation
% This model fit uses data from the first and last trial of each block.
% Reach error is not corrected for baseline error. 

clear all
close all
clc

%% Load the data

data = dlmread('reachError_N4.dat');
data = data(data(:,2) == 1,:);

exList = [118 132];
data(ismember(data(:,1),exList),:) = [];

plist = unique(data(:,1));
nsubs = numel(plist);

disp(aggregate(aggregate(data,[1 3],3,@count),2,2,@count))

%%
[gColor, colorNames] = graphColors(4,0);
tLabels     = {'baseline','low','high','none'};
figTitle    = {'Training 1', 'Training 2', 'Washout', 'Recall' ,'T2 - T1', 'Recall-T2'};
lineWidth   = 2;
markerSize  = 1;
xStart      = [1  41  81 101];
xEnd        = [40 80 100 120]; 
xlimit      = [40 40 20 20];
%%
parms = 0.5;

minParms                    = [-1.0];
maxParms                    = [1.0];
parmSteps                   = [0.10];
psrFactor                   = [0.50];
psaFactor                   = [1.50];

fValue   = nan(nsubs,4);
fParms   = nan(nsubs,4);
rSquared = nan(nsubs,4);

modelPred = nan(nsubs,sum(xlimit));

model = @fit_expcurve;

for pIdx = 1:nsubs
    for cIdx = 1:4
        pData = data(data(:,1) == plist(pIdx) & data(:,4) == cIdx+1,:);
        pno = mean(pData(:,1));
        cnd = mean(pData(:,3));
        
        pData = pData(:,6)';
        nt      = length(pData);
        
        pData(isnan(pData)) = [];
        
        constant = mean(pData(end-5:end));
        delta    = pData(1) - constant;
        
        
        %     options = optimset('Display','final','MaxIter',10000,'MaxFunEvals',1e16,'TolFun',1e-16,'TolX',1e-16);
        %     [fit, fval] = fmincon(model,parms,Aeq,beq,[],[],minParms,maxParms,[],options, [sp delta], pData);
        
        [fit, fval, exitflag] =...
            hkjeeves(model, parms,...
            0, minParms, maxParms,...  % Set first parm to -1 to show iteration
            -1, 1e-4, 10000,...
            parmSteps, psaFactor, psrFactor,...
            [constant delta], pData);
        
        fValue(pIdx,cIdx) = fval;
        fParms(pIdx,cIdx) = fit;
        
        y = expFun([constant delta fit],length(pData));
        
        r = power(corrcoef([pData;y]'),2);
        rSquared(pIdx,cIdx) = r(1,2);
        
        if length(pData) < nt
            y = expFun([constant delta fit],nt);
        end
        
        modelPred(pIdx,xStart(cIdx):xEnd(cIdx)) = y;
    end % for cidx
end % for i...

modelPred = [aggregate(data,[1 3],[1 3],@mean,1), modelPred];

meanPred = aggregate(modelPred,2,3:size(modelPred,2),@nanmean);
meanData = [aggregate(data,3:5,6,@nanmean),aggregate(data,3:5,6,@nanstd,1)./sqrt(10)];

fParms      = [aggregate(data,[1 3],[1 3],@mean,1),abs(fParms)];
fParms(:,7) = fParms(:,4) - fParms(:,3);
fParms(:,8) = fParms(:,6) - fParms(:,4);
meanParms   = aggregate(fParms,2,3:size(fParms,2));
stdParms    = aggregate(fParms,2,3:size(fParms,2),@std,1)./sqrt(10);

%%
figure('windowstyle','docked','color','w')
for i = 1:4
    subplot(1,4,i);
    blkData = meanData(meanData(:,2) == i+1,:);
    
    hold on
    errorbar(blkData(blkData(:,1) == 1,4), blkData(blkData(:,1) == 1,5), ...
        'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
    errorbar(blkData(blkData(:,1) == 2,4), blkData(blkData(:,1) == 2,5), ...
        'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
    errorbar(blkData(blkData(:,1) == 3,4), blkData(blkData(:,1) == 3,5), ...
        'ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
    errorbar(blkData(blkData(:,1) == 4,4), blkData(blkData(:,1) == 4,5), ...
        'ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));

    plot(meanPred(1,(xStart(i):xEnd(i))+1), ...
        '-ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
    plot(meanPred(2,(xStart(i):xEnd(i))+1), ...
        '-ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
    plot(meanPred(3,(xStart(i):xEnd(i))+1), ...
        '-ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
    plot(meanPred(4,(xStart(i):xEnd(i))+1), ...
        '-ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
    hold off

    ylabel('Error (deg)','Fontsize',24)
    xlabel('Trial Blocks','Fontsize',24)
    
    set(gca,'XLim', [.5 xlimit(i)+.5], 'XTick', 0:10:xlimit(i), 'Fontsize', 18)
    set(gca,'YLim', [-40 40],'YTick', -40:10:40);
    title(figTitle{i})
    if i == 1
        legend('baseline','low','high','none');
    end

end

figure('windowstyle','docked','color','w')
for i = 1:6
    subplot(2,3,i);
    hold on
    errorbar(1,meanParms(1,i+1),stdParms(1,i), ...
        'ok','color',gColor(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(1,:));
    errorbar(2,meanParms(2,i+1),stdParms(2,i), ...
        'ok','color',gColor(2,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(2,:));
    errorbar(3,meanParms(3,i+1),stdParms(3,i), ...
        'ok','color',gColor(3,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(3,:));
    errorbar(4,meanParms(4,i+1),stdParms(4,i), ...
        'ok','color',gColor(4,:),'LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',gColor(4,:));
    hold off
    
    ylabel('Learning Rate','Fontsize',24)
    xlabel('Group','Fontsize',24)
    
    set(gca,'XLim', [.5 4.5], 'XTick', 1:4, 'XTickLabel', tLabels, 'Fontsize', 16)
    if i > 4
        set(gca,'YLim', [-.3 0.3],'YTick', -.3:.3:.3,'Fontsize', 16);
    else
        set(gca,'YLim', [0 0.6],'YTick', 0:0.3:0.6,'Fontsize', 16);
    end
    title(figTitle{i})
end % for i...