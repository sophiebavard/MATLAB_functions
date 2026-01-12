function [Nbar,Nsub] = violinplotSB_double(DataCell1,DataCell2,Colors1,Colors2,Yinf,Ysup,varargin)

% Sophie Bavard - July 2024
% Creates 2 violin plots with mean, error bars, confidence interval, kernel density and corresponding datapoints.

% Parse optional inputs
p = inputParser;
addOptional(p, 'DotSize', 37.5, @isnumeric);  % default is Wbar*50 where Wbar=0.75
parse(p, varargin{:});

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
% dot size - from input parameter
Wdot = p.Results.DotSize;
% middle space
space = Wbar/2.5;

% confidence interval
ConfInter = 0.95;

% Get current figure handle
fig = gcf;
ax = gca;

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
    if abs(max(DataMatrix1) - min(DataMatrix1)) < 1e-10  % all values identical
        value1   = 1;      % just one vertical position
        density1 = 1;      % fake flat density for plotting
    else
        if iqr(DataMatrix1) ~= 0
            [density1, value1] = ksdensity(DataMatrix1, 'Bandwidth', 0.9 * min(std(DataMatrix1), iqr(DataMatrix1)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
        else
            [density1, value1] = ksdensity(DataMatrix1, 'Bandwidth', 0.9 * std(DataMatrix1) * Nsub^(-1/5));
        end
        density1 = density1(value1 >= min(DataMatrix1) & value1 <= max(DataMatrix1));
        value1 = value1(value1 >= min(DataMatrix1) & value1 <= max(DataMatrix1));
        value1(1) = min(DataMatrix1);
        value1(end) = max(DataMatrix1);
    end

    width1 = Wbar/2/max(density1);

    % calculate kernel density estimation for the violin
    if abs(max(DataMatrix2) - min(DataMatrix2)) < 1e-10  % all values identical
        value2   = 1;      % just one vertical position
        density2 = 1;      % fake flat density for plotting
    else
        if iqr(DataMatrix2) ~= 0
            [density2, value2] = ksdensity(DataMatrix2, 'Bandwidth', 0.9 * min(std(DataMatrix2), iqr(DataMatrix2)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
        else
            [density2, value2] = ksdensity(DataMatrix2, 'Bandwidth', 0.9 * std(DataMatrix2) * Nsub^(-1/5));
        end
        density2 = density2(value2 >= min(DataMatrix2) & value2 <= max(DataMatrix2));
        value2 = value2(value2 >= min(DataMatrix2) & value2 <= max(DataMatrix2));
        value2(1) = min(DataMatrix2);
        value2(end) = max(DataMatrix2);
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
        d = interp1(value1, density1*width1, [curve1-sem1*conf:step:curve1+sem1*conf curve1+sem1*conf]);
    fill([n n-d n] - space,...
        [curve1-sem1*conf curve1-sem1*conf:step:curve1+sem1*conf curve1+sem1*conf curve1+sem1*conf],...
        Colors1(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.25);
    end
    hold on
    if length(density2) > 1
        d = interp1(value2, density2*width2, [curve2-sem2*conf:step:curve2+sem2*conf curve2+sem2*conf]);
    fill([n n+d n] + space,...
        [curve2-sem2*conf curve2-sem2*conf:step:curve2+sem2*conf curve2+sem2*conf curve2+sem2*conf],...
        Colors2(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.25);
    end
    hold on
    
    % ERROR BAR INTERVALS   
    if length(density1) > 1
        d = interp1(value1, density1*width1, [curve1-sem1:step:curve1+sem1 curve1+sem1]);
    fill([n n-d n] - space,...
        [curve1-sem1 curve1-sem1:step:curve1+sem1 curve1+sem1 curve1+sem1],...
        Colors1(n,:),...
        'EdgeColor', 'none',...
        'FaceAlpha',0.6);
    end
    hold on   
    if length(density2) > 1
        d = interp1(value2, density2*width2, [curve2-sem2:step:curve2+sem2 curve2+sem2]);
    fill([n n+d n] + space,...
        [curve2-sem2 curve2-sem2:step:curve2+sem2 curve2+sem2 curve2+sem2],...
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
    scatter(n + jitter1 - space + space/5, DataMatrix1, Wdot,...
        Colors1(n,:),'filled',...
        'marker','o',...
        'MarkerFaceAlpha',0.4);
    hold on
    jitter2=(Wbar/8)*abs(zscore(1:length(DataMatrix2))'/max(zscore(1:length(DataMatrix2))'));
    scatter(n - jitter2 + space - space/5, DataMatrix2, Wdot,...
        Colors2(n,:),'filled',...
        'marker','o',...
        'MarkerFaceAlpha',0.4);
    hold on

    % LINES BETWEEN DOTS -- GREY
    % plot([n + jitter1 - space + space/8 n - jitter2 + space - space/8]',...
    %     [DataMatrix1 DataMatrix2]',...
    %     'Color',[0 0 0 0.2]);
    % hold on

    % LINES BETWEEN DOTS -- COLOR GRADIENT
    for i = 1:length(DataMatrix1)
        patch('XData',[n + jitter1(i) - space + space/5, n - jitter2(i) + space - space/5], ...
            'YData',[DataMatrix1(i), DataMatrix2(i)], ...
            'ZData',[0 0], ...                         % dummy z
            'FaceColor','none', ...                    % no fill
            'EdgeColor','interp', ...                  % interpolate along edge
            'EdgeAlpha',0.2, ...                       % line transparency
            'LineWidth',0.5, ...                       % line size
            'FaceVertexCData',[Colors1(n,:);Colors2(n,:)]);      % colors per vertex
    end
    hold on

end

% axes and stuff

set(gca,'XTickLabel',[])
set(gca,'XTick',[])

ylim([Yinf Ysup]);
xlim([0+Wbar/2 Nbar+1-Wbar/2]);

% DYNAMIC DOT SIZE - Store reference and setup callback

% Get all scatter objects and store their original sizes
scatterObjs = findobj(ax, 'Type', 'scatter');
for i = 1:length(scatterObjs)
    userData = get(scatterObjs(i), 'UserData');
    userData.originalSize = get(scatterObjs(i), 'SizeData');
    set(scatterObjs(i), 'UserData', userData);
end

% Store initial axes size in pixels for reference
axPos = getpixelposition(ax);
axUserData = get(ax, 'UserData');
axUserData.referenceSize = sqrt(axPos(3) * axPos(4));
axUserData.hasViolinPlots = true;
set(ax, 'UserData', axUserData);

% Get or initialize the list of violin axes in the figure
figUserData = get(fig, 'UserData');
if ~isfield(figUserData, 'violinAxes') || isempty(figUserData.violinAxes)
    figUserData.violinAxes = ax;
else
    % Check if this axes is already in the list
    alreadyExists = false;
    for i = 1:length(figUserData.violinAxes)
        if isequal(figUserData.violinAxes(i), ax)
            alreadyExists = true;
            break;
        end
    end
    % Add this axes to the list if not already present
    if ~alreadyExists
        figUserData.violinAxes(end+1) = ax;
    end
end
set(fig, 'UserData', figUserData);

% Set up dynamic resizing callback (only once per figure)
currentCallback = get(fig, 'SizeChangedFcn');
if isempty(currentCallback)
    set(fig, 'SizeChangedFcn', @(src,evt) updateAllViolinAxes(src));
end

end


function updateAllViolinAxes(fig)
% Update all axes in the figure that have violin plots
try
    figUserData = get(fig, 'UserData');
    if ~isfield(figUserData, 'violinAxes')
        return;
    end
    
    allAxes = figUserData.violinAxes;
    
    % Create list of valid axes
    validAxes = [];
    for i = 1:length(allAxes)
        if ishghandle(allAxes(i), 'axes')
            validAxes(end+1) = allAxes(i);
        end
    end
    
    % Update the list
    figUserData.violinAxes = validAxes;
    set(fig, 'UserData', figUserData);
    
    % Update each valid axes
    for i = 1:length(validAxes)
        updateMarkerSizes(validAxes(i));
    end
catch ME
    warning('Error in updateAllViolinAxes: %s', ME.message);
end
end


function updateMarkerSizes(ax)
% Get current and reference axes sizes
% Check if axes is still valid
if ~ishghandle(ax, 'axes')
    return;
end

axUserData = get(ax, 'UserData');
if ~isfield(axUserData, 'referenceSize')
    return;
end

axPos = getpixelposition(ax);
currentSize = sqrt(axPos(3) * axPos(4));
referenceSize = axUserData.referenceSize;

% Calculate scaling factor
scaleFactor = currentSize / referenceSize;

% Find all scatter objects in the axes
scatterObjs = findobj(ax, 'Type', 'scatter');

% Update all scatter objects
for i = 1:length(scatterObjs)
    if ishghandle(scatterObjs(i))
        scatterUserData = get(scatterObjs(i), 'UserData');
        if isfield(scatterUserData, 'originalSize')
            originalSize = scatterUserData.originalSize;
            set(scatterObjs(i), 'SizeData', originalSize * scaleFactor);
        end
    end
end
end