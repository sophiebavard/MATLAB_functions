function [Nbar,Nsub] = violinplotSB(DataCell,Colors,Yinf,Ysup,varargin)

% Sophie Bavard - January 2024
% Creates a violin plot with mean, error bars, confidence interval, kernel density.

% Parse optional inputs
p = inputParser;
addOptional(p, 'DotSize', 37.5, @isnumeric);  % default is Wbar*50 where Wbar=0.75
parse(p, varargin{:});

% transforms the Data matrix into cell format if needed
if iscell(DataCell)==0
    DataCell = num2cell(DataCell,2);
end

% number of factors/groups/conditions
Nbar = size(DataCell,1);
% bar size
Wbar = 0.75;
% dot size - from input parameter
Wdot = p.Results.DotSize;

% confidence interval
ConfInter = 0.95;

% Get current figure handle
fig = gcf;
ax = gca;

for n = 1:Nbar

    clear DataMatrix
    clear jitter jitterstrength
    DataMatrix = DataCell{n,:}';

    % Find the step to adapt the plots (maximum nb of digits)
    scale = ceil(log10(max(abs(DataMatrix))));
    step  = 10^scale/1000;

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
    if abs(max(DataMatrix) - min(DataMatrix)) < 1e-10  % all values identical
        value   = 1;      % just one vertical position
        density = 1;      % fake flat density for plotting
    else
        if iqr(DataMatrix) ~= 0
            [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * min(std(DataMatrix), iqr(DataMatrix)/1.34) * Nsub^(-1/5)); % change Bandwidth for violin shape. Default MATLAB: std(DataMatrix)*(4/(3*Nsub))^(1/5)
        else
            [density, value] = ksdensity(DataMatrix, 'Bandwidth', 0.9 * std(DataMatrix) * Nsub^(-1/5));
        end
        density = density(value >= min(DataMatrix) & value <= max(DataMatrix));
        value = value(value >= min(DataMatrix) & value <= max(DataMatrix));
        value(1) = min(DataMatrix);
        value(end) = max(DataMatrix);
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
        d = interp1(value, density*width, [curve-sem*conf:step:curve+sem*conf curve+sem*conf]);
        fill([n n+d n],...
            [curve-sem*conf curve-sem*conf:step:curve+sem*conf curve+sem*conf curve+sem*conf],...
            Colors(n,:),...
            'EdgeColor', 'none',...
            'FaceAlpha',0.25);
    end
    hold on

    % ERROR BAR INTERVAL
    if length(density) > 1
        d = interp1(value, density*width, [curve-sem:step:curve+sem curve+sem]);
        fill([n n+d n],...
            [curve-sem curve-sem:step:curve+sem curve+sem curve+sem],...
            Colors(n,:),...
            'EdgeColor', 'none',...
            'FaceAlpha',0.6);
    end
    hold on

    % INDIVIDUAL DOTS with dynamic scaling
    jitter = abs(zscore(1:length(DataMatrix))'/max(zscore(1:length(DataMatrix))'));

    % Create scatter plot with transparency
    scatter(n - Wbar/10 - jitter.*(Wbar/2- Wbar/10), DataMatrix, Wdot,...
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
    plot(xMean,yMean,'-','LineWidth',1.5,'Color','k');
    hold on

end

% axes and stuff
set(gca,'XTickLabel',[])
set(gca,'XTick',[])
ylim([Yinf Ysup]);
xlim([0 Nbar+1]);

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