function [Nbar,Nsub] = violinplotSB(DataCell,Colors,Yinf,Ysup)

% Sophie Bavard - January 2024
% Creates a violin plot with mean, error bars, confidence interval, kernel density.

% transforms the Data matrix into cell format if needed
if iscell(DataCell)==0
    DataCell = num2cell(DataCell,2);
end

% number of factors/groups/conditions
Nbar = size(DataCell,1);
% bar size
Wbar = 0.75;

% confidence interval
ConfInter = 0.95;

for n = 1:Nbar
    
    clear DataMatrix
    clear jitter jitterstrength
    DataMatrix = DataCell{n,:}';

    % if all NaNs
    if sum(isnan(DataMatrix))==size(DataMatrix,1)
        DataMatrix = 0;
    end
    
    % number of subjects
    Nsub = length(DataMatrix(~isnan(DataMatrix)));
    
    curve = nanmean(DataMatrix);
    sem   = nanstd(DataMatrix')'/sqrt(Nsub);
    conf  = tinv(1 - 0.5*(1-ConfInter),Nsub);
    
    % PLOT THE VIOLINS
    
    % calculate kernel density estimation for the violin
    if iqr(DataMatrix) ~= 0
        [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * min(std(DataMatrix), iqr(DataMatrix)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
    else
        [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * std(DataMatrix) * Nsub^(-1/5));
    end
    density = density(value >= min(DataMatrix) & value <= max(DataMatrix));
    value = value(value >= min(DataMatrix) & value <= max(DataMatrix));
    value(1) = min(DataMatrix);
    value(end) = max(DataMatrix);
    
    % all data is identical
    if min(DataMatrix) == max(DataMatrix)
        density = 1; value = 1;
    end
    width = Wbar/2/max(density);
    
    % plot the violin
    fill([n n+density*width n],...
        [value(1) value value(end)],...
        Colors(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.2);
    hold on
    
    % CONFIDENCE INTERVAL   
    if length(density) > 1
        d = interp1(value, density*width, curve-sem*conf:0.0001:curve+sem*conf);
    fill([n n+d n],...
        [curve-sem*conf curve-sem*conf:0.0001:curve+sem*conf curve+sem*conf],...
        Colors(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.25);
    end
    hold on
    
    % ERROR BAR INTERVAL   
    if length(density) > 1
        d = interp1(value, density*width, curve-sem:0.0001:curve+sem);
    fill([n n+d n],...
        [curve-sem curve-sem:0.0001:curve+sem curve+sem],...
        Colors(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.6);
    end
    hold on
        
    % INDIVIDUAL DOTS
    if length(density) > 1
        jitterstrength = interp1(value, density*width, DataMatrix);
    else % all data is identical
        jitterstrength = density*width;
    end

    jitter=abs(zscore(1:length(DataMatrix))'/max(zscore(1:length(DataMatrix))'));
    
	scatter(n - Wbar/10 - jitter.*(Wbar/2- Wbar/10), DataMatrix, 10,...
        Colors(n,:),'filled',...
        'marker','o',...
        'MarkerFaceAlpha',0.4);
    hold on
    
    % MEAN HORIZONTAL BAR
    if length(density)>1
        xM = interp1(value, density*width, curve);
    else
        xM = density*width;
    end
    xMean = [n ; n + xM];
    yMean = [curve; curve];
    plot(xMean,yMean,'-','LineWidth',1,'Color','k');
    hold on
    
end

% axes and stuff

set(gca,'XTickLabel',[])
set(gca,'XTick',[])

ylim([Yinf Ysup]);
xlim([0 Nbar+1]);












