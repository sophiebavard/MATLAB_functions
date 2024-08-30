function [Nbar,Nsub] = violinplotSB_double(DataCell1,DataCell2,Colors1,Colors2,Yinf,Ysup)

% Sophie Bavard - July 2024
% Creates 2 violin plots with mean, error bars, confidence interval, kernel density and corresponding datapoints.

% transforms the Data matrix into cell format if needed
if iscell(DataCell1)==0
    DataCell1 = num2cell(DataCell1,2);
end
if iscell(DataCell2)==0
    DataCell2 = num2cell(DataCell2,2);
end

% number of factors/groups/conditions
Nbar = size(DataCell1,1);
% bar size
Wbar = 0.5;
% dot size
Wdot = 10;
% middle space
space = Wbar/2.5;

% confidence interval
ConfInter = 0.95;

for n = 1:Nbar
    
    clear DataMatrix1
    clear DataMatrix2
  
    DataMatrix1 = DataCell1{n,:}';
    DataMatrix2 = DataCell2{n,:}';

    % Find the step to adapt the plots (maximum nb of digits)
    scale = ceil(log10(max(max(abs(DataMatrix1),abs(DataMatrix2)))));
    step  = 10^scale/1000;

    % if all NaNs
    if sum(isnan(DataMatrix1))==size(DataMatrix1,1)
        DataMatrix1 = 0;
    end
    
    % number of subjects
    Nsub = length(DataMatrix1(~isnan(DataMatrix1)));
    conf  = tinv(1 - 0.5*(1-ConfInter),Nsub);

    curve1 = nanmean(DataMatrix1);
    sem1   = nanstd(DataMatrix1')'/sqrt(Nsub);

    curve2 = nanmean(DataMatrix2);
    sem2   = nanstd(DataMatrix2')'/sqrt(Nsub);

    % PLOT THE VIOLINS
    
    % calculate kernel density estimation for the violin
    if iqr(DataMatrix1) ~= 0
        [density1, value1] = ksdensity(DataMatrix1, 'Bandwidth', 0.9 * min(std(DataMatrix1), iqr(DataMatrix1)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
    else
        [density1, value1] = ksdensity(DataMatrix1, 'Bandwidth', 0.9 * std(DataMatrix1) * Nsub^(-1/5));
    end
    density1 = density1(value1 >= min(DataMatrix1) & value1 <= max(DataMatrix1));
    value1 = value1(value1 >= min(DataMatrix1) & value1 <= max(DataMatrix1));
    value1(1) = min(DataMatrix1);
    value1(end) = max(DataMatrix1);
    % all data is identical
    if min(DataMatrix1) == max(DataMatrix1)
        density1 = 1; value1 = 1;
    end
    width1 = Wbar/2/max(density1);
    
    % calculate kernel density estimation for the violin
    if iqr(DataMatrix2) ~= 0
        [density2, value2] = ksdensity(DataMatrix2, 'Bandwidth', 0.9 * min(std(DataMatrix2), iqr(DataMatrix2)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
    else
        [density2, value2] = ksdensity(DataMatrix2, 'Bandwidth', 0.9 * std(DataMatrix2) * Nsub^(-1/5));
    end
    density2 = density2(value2 >= min(DataMatrix2) & value2 <= max(DataMatrix2));
    value2 = value2(value2 >= min(DataMatrix2) & value2 <= max(DataMatrix2));
    value2(1) = min(DataMatrix2);
    value2(end) = max(DataMatrix2);
    % all data is identical
    if min(DataMatrix2) == max(DataMatrix2)
        density2 = 1; value2 = 1;
    end
    width2 = Wbar/2/max(density2);
    
    % VIOLINS
    fill([n n-density1*width1 n] - space,...
        [value1(1) value1 value1(end)],...
        Colors1(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.2);
    hold on
    fill([n n+density2*width2 n] + space,...
        [value2(1) value2 value2(end)],...
        Colors2(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.2);
    hold on
    
    % CONFIDENCE INTERVALS
    if length(density1) > 1
        d = interp1(value1, density1*width1, curve1-sem1*conf:step:curve1+sem1*conf);
    fill([n n-d n] - space,...
        [curve1-sem1*conf curve1-sem1*conf:step:curve1+sem1*conf curve1+sem1*conf],...
        Colors1(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.25);
    end
    hold on
    if length(density2) > 1
        d = interp1(value2, density2*width2, curve2-sem2*conf:step:curve2+sem2*conf);
    fill([n n+d n] + space,...
        [curve2-sem2*conf curve2-sem2*conf:step:curve2+sem2*conf curve2+sem2*conf],...
        Colors2(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.25);
    end
    hold on
    
    % ERROR BAR INTERVALS   
    if length(density1) > 1
        d = interp1(value1, density1*width1, curve1-sem1:step:curve1+sem1);
    fill([n n-d n] - space,...
        [curve1-sem1 curve1-sem1:step:curve1+sem1 curve1+sem1],...
        Colors1(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.6);
    end
    hold on   
    if length(density2) > 1
        d = interp1(value2, density2*width2, curve2-sem2:step:curve2+sem2);
    fill([n n+d n] + space,...
        [curve2-sem2 curve2-sem2:step:curve2+sem2 curve2+sem2],...
        Colors2(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.6);
    end
    hold on
    
    % MEAN HORIZONTAL BARS
    if length(density1)>1
        xM1 = interp1(value1, density1*width1, curve1);
    else
        xM1 = density1*width1;
    end
    xMean1 = [n ; n - xM1] - space;
    yMean1 = [curve1; curve1];
    plot(xMean1,yMean1,'-','LineWidth',1.5,'Color','k');
    hold on
    if length(density2)>1
        xM2 = interp1(value2, density2*width2, curve2);
    else
        xM2 = density2*width2;
    end
    xMean2 = [n ; n + xM2] + space;
    yMean2 = [curve2; curve2];
    plot(xMean2,yMean2,'-','LineWidth',1.5,'Color','k');
    hold on
        
    % INDIVIDUAL DOTS
    jitter1=(Wbar/8)*abs(zscore(1:length(DataMatrix1))'/max(zscore(1:length(DataMatrix1))'));
    scatter(n + jitter1 - space + space/8, DataMatrix1, Wdot,...
        Colors1(n,:),'filled',...
        'marker','o',...
        'MarkerFaceAlpha',0.4);
    hold on
    jitter2=(Wbar/8)*abs(zscore(1:length(DataMatrix2))'/max(zscore(1:length(DataMatrix2))'));
    scatter(n - jitter2 + space - space/8, DataMatrix2, Wdot,...
        Colors2(n,:),'filled',...
        'marker','o',...
        'MarkerFaceAlpha',0.4);
    hold on

    % LINES BETWEEN DOTS
    plot([n + jitter1 - space + space/8 n - jitter2 + space - space/8]',...
        [DataMatrix1 DataMatrix2]',...
        'Color',[0 0 0 0.2]);
    hold on
    
end

% axes and stuff

set(gca,'XTickLabel',[])
set(gca,'XTick',[])

ylim([Yinf Ysup]);
xlim([0+Wbar/2 Nbar+1-Wbar/2]);

% Check potential plotting issues
% profile on
% figure; ...
% profile viewer












