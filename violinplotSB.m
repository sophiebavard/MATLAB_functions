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
% dot size
Wdot = Wbar*50;

% confidence interval
ConfInter = 0.95;

% Get current figure handle
fig = gcf;

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
    plot(xMean,yMean,'-','LineWidth',1,'Color','k');
    hold on

end

% axes and stuff
set(gca,'XTickLabel',[])
set(gca,'XTick',[])
ylim([Yinf Ysup]);
xlim([0 Nbar+1]);

% DYNAMIC DOT SIZE

% Set up dynamic resizing callback
fig.SizeChangedFcn = @(src,evt) updateAllMarkerSizes(fig);

% Store initial axes size for reference scaling
ax = gca;
axPos = getpixelposition(ax);
referenceSize = sqrt(axPos(3) * axPos(4));

% Store reference size in axes UserData for scaling calculations
ax.UserData.referenceSize = referenceSize;

end


function updateAllMarkerSizes(fig)
% Find all axes in the figure
allAxes = findobj(fig, 'Type', 'axes');

% Update marker sizes for each axes
for i = 1:length(allAxes)
    updateMarkerSizes(allAxes(i));
end
end

function updateMarkerSizes(ax)
% Get current axes position in pixels
try
    axPos = getpixelposition(ax);
    figPos = getpixelposition(get(ax, 'Parent'));

    % Get normalized reference size (initial relative size when plot was created)
    if isfield(ax.UserData, 'normalizedReferenceSize')
        normalizedReferenceSize = ax.UserData.normalizedReferenceSize;
    else
        % If no reference stored, don't scale
        return;
    end

    % Calculate current normalized size (axes size relative to current figure size)
    axesArea = axPos(3) * axPos(4);
    figureArea = figPos(3) * figPos(4);
    currentNormalizedSize = sqrt(axesArea / figureArea);

    % Calculate scaling factor relative to normalized sizes
    scaleFactor = currentNormalizedSize / normalizedReferenceSize;

    % Find all scatter objects in the axes
    scatterObjs = findobj(ax, 'Type', 'scatter');

    % Update all scatter objects
    for i = 1:length(scatterObjs)
        if isvalid(scatterObjs(i))
            % Get original size from UserData, or store it if first time
            if isempty(scatterObjs(i).UserData)
                % Store the original size (Wbar*25) on first scaling call
                scatterObjs(i).UserData = scatterObjs(i).SizeData;
            end
            originalSize = scatterObjs(i).UserData;

            % Apply scaling while maintaining original size as baseline
            scatterObjs(i).SizeData = originalSize * scaleFactor;
        end
    end
catch
    % Handle any errors silently
end
end




