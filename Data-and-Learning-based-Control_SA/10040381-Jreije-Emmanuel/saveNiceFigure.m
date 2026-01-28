function saveNiceFigure(filename, layoutType)
% saveNiceFigure Formats the current figure for A4 and saves it.
%
%   saveNiceFigure(filename)
%   saveNiceFigure(filename, layoutType)
%
%   filename   : base name without extension (e.g. 'fig_inputs')
%   layoutType : 'single', 'wide', or 'tall'  (optional, default: 'wide')

    if nargin < 2
        layoutType = 'wide';  % default for your 3-subplot figure
    end

    fig = gcf;  % current figure handle

    % ---- layout options in centimeters ----
    switch lower(layoutType)
        case 'single'
            width  = 15;  % cm
            height = 9;   % cm
        case 'wide'
            width  = 22;  % cm
            height = 14;  % cm
        case 'tall'
            width  = 15;  % cm
            height = 18;  % cm
        otherwise
            warning('Unknown layoutType. Using ''wide''.');
            width  = 18;
            height = 14;
    end

    % ---- set figure size ----
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [2 2 width height]);  % [left bottom width height]

    % ---- global font settings (you can tune these) ----
    font_title = 14;
    font_label = 14;
    font_legend = 10;
    font_ticks = 10;

    % apply to all axes
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k), 'FontSize', font_ticks);               % tick labels
        % if labels already exist, update them
        set(get(ax(k),'XLabel'), 'FontSize', font_label);
        set(get(ax(k),'YLabel'), 'FontSize', font_label);
        set(get(ax(k),'Title'),  'FontSize', font_title);
    end

    % apply to all legends
    lg = findall(fig, 'Type', 'legend');
    for k = 1:numel(lg)
        set(lg(k), 'FontSize', font_legend);
    end

    % ---- Paper settings for export ----
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [width height]);
    set(fig, 'PaperPosition', [0 0 width height]);
    set(fig, 'PaperPositionMode', 'manual');

    % ---- export ----
    dpi = 300;

    % PDF
    print(fig, [filename '.pdf'], '-dpdf', sprintf('-r%d', dpi));

    % PNG
    print(fig, [filename '.png'], '-dpng', sprintf('-r%d', dpi));

    fprintf('Figure saved as %s.pdf and %s.png\n', filename, filename);
end
