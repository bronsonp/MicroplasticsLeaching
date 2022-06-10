% This code makes the plots that are used in the paper. First run the
% `curve_fitting` script to produce the results, then run this script to
% produce the figures.

%% Load data 
load("Results\data.mat");
fits = struct();
for p = ["BPS", "BPA", "DEHT", "DEHP"]
    fits.(p) = load(fullfile("Results", p+".mat"));
end

colours = {...
    [0.21484375,0.4921875,0.71875], ...
    [0.30078125,0.68359375,0.2890625], ...
    [0.890625,0.1015625,0.109375], ...
    [0.59375,0.3046875,0.63671875], ...
    [0.99609375,0.49609375,0], ...
    [0.6484375,0.3359375,0.15625], ...
    [0.96484375,0.50390625,0.74609375], ...
    [0.5,0.5,0.5], ...
    [0.8965,0.8965,0.1793]};
markers = {'o', '^', 's', 'd'};
font = 'Segoe UI';


%% Table of values
f = fits.BPA;
for p = ["BPS", "DEHT", "DEHP"]
    f(end+1) = fits.(p);
end
f = struct2table(f);
% convert units
f.D_plastic_25C = f.D_plastic_25C * (10^-4)^2 / 86400; % cm^2/s
f.D_plastic_EA = f.D_plastic_EA / 1000; % kJ/mol
f.delta_still = f.delta_still * 1e-6 * 1e3; % mm
f.delta_agitated = f.delta_agitated * 1e-6 * 1e3; % mm
f






%% Figure showing varying temperature

figure();
sz = "<200um";
wc = "Still";
plot_data = data(data.Size == sz & data.WaterCondition == wc, :);
subplot_id = 1;
tiles=tiledlayout(2,2,'TileSpacing', 'compact');
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    nexttile();
    line_id = 1;
    lg = [];
    
    for temp = unique(plot_data.Temperature)'
        this_data = plot_data(plot_data.Plasticiser == p & plot_data.Temperature == temp, :);
        if ~isempty(this_data)
            % plot the data, including the pre-rinsed amount
            time = [-4; this_data.Time];
            wtpc = [this_data.MeanPlasticiserWtPercent(1); this_data.InternalPlasticiserWtPercent];
            wtpc_std = [this_data.StdPlasticiserWtPercent(1); this_data.StdPlasticiserWtPercent];
            
            lg(end+1) = errorbar(time, wtpc, wtpc_std, markers{line_id}, ...
                'Color', colours{line_id}, ...
                'MarkerEdgeColor', colours{line_id}, ...
                'MarkerFaceColor', colours{line_id}, ...
                'DisplayName', sprintf('%g °C', temp));
            hold on;
        
            % plot the model
            p0 = this_data.InternalPlasticiserConc(this_data.Time == 0, :);
            [time, model_wtpercent] = predict(fits.(p), this_data);
            plot(time, model_wtpercent, 'Color', colours{line_id});
            
            line_id = line_id + 1;
            
        end
    end

    ax = gca();
    ax.FontName = font;
    ax.FontSize = 11;
    ax.XLim = [-4 22];
    ax.XTick = [-4 0:5:20];
    ax.XTickLabel{1} = 'BR';
    ax.TickLength = [0.02 0.05];
%     ax.TickDir = 'out';
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";

% add broken axis annotation
annotations = [];
line_len = 0.03;
line_sep = 0.01;
left_offset = 0.01;
rot = 20;

for child = tiles.Children'
    if class(child) ~= "matlab.graphics.axis.Axes"
        continue;
    end

    left = child.Position(1)+left_offset;
    for bottom = [child.Position(2)-line_len/2, child.Position(2)+child.Position(4)-line_len/2]
        annotations(end+1) = annotation("rectangle", [left, bottom, line_sep+0.0001, line_len], "FaceColor", "w", "Rotation", -rot, "Color", "none");
        annotations(end+1) = annotation("line", left + [0 line_len*sind(rot)], bottom + [0 line_len*cosd(rot)]);
        annotations(end+1) = annotation("line", left + line_sep + [0 line_len*sind(rot)], bottom + [0 line_len*cosd(rot)]);
    end
end


%% Figure showing varying agitation

figure();
sz = "<200um";
temp = 25.5;
plot_data = data(data.Size == sz & data.Temperature == temp, :);
subplot_id = 1;
tiles=tiledlayout(2,2,'TileSpacing', 'compact');
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    nexttile();
    line_id = 1;
    lg = [];
    
    for wc = unique(plot_data.WaterCondition)'
        this_data = plot_data(plot_data.Plasticiser == p & plot_data.WaterCondition == wc, :);
        if ~isempty(this_data)
            % plot the data, including the pre-rinsed amount
            time = [-4; this_data.Time];
            wtpc = [this_data.MeanPlasticiserWtPercent(1); this_data.InternalPlasticiserWtPercent];
            wtpc_std = [this_data.StdPlasticiserWtPercent(1); this_data.StdPlasticiserWtPercent];
            
            if wc == "Still"
                mfc = "none";
            else
                mfc = colours{line_id};
            end
            lg(end+1) = errorbar(time, wtpc, wtpc_std, markers{line_id}, ...
                'Color', colours{line_id}, ...
                'MarkerEdgeColor', colours{line_id}, ...
                'MarkerFaceColor', mfc, ...
                'DisplayName', sprintf('%s', wc));
            hold on;
        
            % plot the model
            p0 = this_data.InternalPlasticiserConc(this_data.Time == 0, :);
            [time, model_wtpercent] = predict(fits.(p), this_data);
            plot(time, model_wtpercent, 'Color', colours{line_id});
            
            line_id = line_id + 1;
            
        end
    end

    ax = gca();
    ax.FontName = font;
    ax.FontSize = 11;
    ax.XLim = [-4 22];
    ax.XTick = [-4 0:5:20];
    ax.XTickLabel{1} = 'BR';
    ax.TickLength = [0.02 0.05];
%     ax.TickDir = 'out';
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";

% add broken axis annotation
annotations = [];
line_len = 0.03;
line_sep = 0.01;
left_offset = 0.01;
rot = 20;

for child = tiles.Children'
    if class(child) ~= "matlab.graphics.axis.Axes"
        continue;
    end

    left = child.Position(1)+left_offset;
    for bottom = [child.Position(2)-line_len/2, child.Position(2)+child.Position(4)-line_len/2]
        annotations(end+1) = annotation("rectangle", [left, bottom, line_sep+0.0001, line_len], "FaceColor", "w", "Rotation", -rot, "Color", "none");
        annotations(end+1) = annotation("line", left + [0 line_len*sind(rot)], bottom + [0 line_len*cosd(rot)]);
        annotations(end+1) = annotation("line", left + line_sep + [0 line_len*sind(rot)], bottom + [0 line_len*cosd(rot)]);
    end
end

%% Evaluate the spatial behaviour of each plasticiser

figure();
sz = "<200um";
temp = 25.5;
wc = "Agitated";
plot_data = data(data.Size == sz & data.Temperature == temp & data.WaterCondition == wc, :);
subplot_id = 1;
tiles=tiledlayout(2,2,'TileSpacing', 'compact');

for p = ["BPA", "BPS", "DEHT", "DEHP"]
    this_data = plot_data(plot_data.Plasticiser == p & plot_data.WaterCondition == wc, :);

    % evaluate the model and find the time it takes to decay to 1% of the
    % starting value
    time = linspace(0, 60, 1000);
    [~, ~, y, associated] = predict(fits.(p), this_data, time);
    y = y - associated;
    y_scale = y(1);
    y = y / y_scale;

    target = 0.01;
    mask = y > 1e-4;
    time = time(mask);
    y = y(mask);
    t_leach = interp1(y, time, 0.01)
    times = t_leach .* [0 0.01 0.05 0.1:0.1:1];

    [~, ~, ~, ~, free_plasticiser, r] = predict(fits.(p), this_data, times);
    free_plasticiser = free_plasticiser ./ y_scale;

    
   % make the plot
    nexttile
    lg = [];
    for i = 1:numel(times)
        hsv = [(i-1)/(numel(times)) 0.9 0.8];
        rgb = hsv2rgb(hsv);
        lg(end+1) = plot(r, free_plasticiser(:, i), 'Color', rgb, 'LineWidth', 1, 'DisplayName', sprintf('t / t_{leach} = %g', times(i)/t_leach));
        hold on;
    end
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end
xlabel(tiles, 'Radial position (μm)');
ylabel(tiles, 'Relative concentration of plasticiser, p / p_0');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";


%% Figure showing varying sizes

figure();
temp = 25.5;
wc = "Still";
plot_data = data(data.Temperature == temp & data.WaterCondition == wc, :);
subplot_id = 1;
tiles=tiledlayout(2,1,'TileSpacing', 'compact');
for p = ["BPA", "DEHT"]
    nexttile();
    line_id = 1;
    lg = [];
    
    for sz = unique(plot_data.Size)'
        this_data = plot_data(plot_data.Plasticiser == p & plot_data.Size == sz, :);
        if ~isempty(this_data)
            % plot the data
            time = [this_data.Time];
            wtpc = [this_data.InternalPlasticiserWtPercent];
            wtpc_std = [this_data.StdPlasticiserWtPercent];
            
            lg(end+1) = errorbar(time, wtpc, wtpc_std, markers{line_id}, ...
                'Color', colours{line_id}, ...
                'MarkerEdgeColor', colours{line_id}, ...
                'MarkerFaceColor', colours{line_id}, ...
                'DisplayName', sprintf('%s', sz));
            hold on;
        
            % plot the model
            p0 = this_data.InternalPlasticiserConc(this_data.Time == 0, :);
            [time, model_wtpercent] = predict(fits.(p), this_data);
            plot(time, model_wtpercent, 'Color', colours{line_id});
            
            line_id = line_id + 1;
            
        end
    end

    ax = gca();
    ax.FontName = font;
    ax.FontSize = 11;
    ax.XLim = [0 22];
    ax.XTick = [0:5:20];
    ax.TickLength = [0.02 0.05];
%     ax.TickDir = 'out';
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(1), lg);
l.Location = "northeastoutside";


%% What boundary layer thickness would limit in still water?

figure();
sz = "<200um";
temp = 25.5;
plot_data = data(data.Size == sz & data.Temperature == temp, :);
subplot_id = 1;
tiles=tiledlayout('flow', 'TileSpacing', 'compact');

adjusted_fits = fits;
adjusted_fits.DEHT.delta_agitated = 100 * 1e3; 
adjusted_fits.DEHP.delta_agitated = 100 * 1e3; 


for p = ["DEHT", "DEHP"]
    nexttile();
    line_id = 1;
    lg = [];
    
    for wc = unique(plot_data.WaterCondition)'
        this_data = plot_data(plot_data.Plasticiser == p & plot_data.WaterCondition == wc, :);
        if ~isempty(this_data)
            % plot the data, including the pre-rinsed amount
            time = [-4; this_data.Time];
            wtpc = [this_data.MeanPlasticiserWtPercent(1); this_data.InternalPlasticiserWtPercent];
            wtpc_std = [this_data.StdPlasticiserWtPercent(1); this_data.StdPlasticiserWtPercent];
            
            if wc == "Still"
                mfc = "none";
            else
                mfc = colours{line_id};
            end
            lg(end+1) = errorbar(time, wtpc, wtpc_std, markers{line_id}, ...
                'Color', colours{line_id}, ...
                'MarkerEdgeColor', colours{line_id}, ...
                'MarkerFaceColor', mfc, ...
                'DisplayName', sprintf('%s', wc));
            hold on;
        
            % plot the model
            p0 = this_data.InternalPlasticiserConc(this_data.Time == 0, :);


            [time, model_wtpercent] = predict(adjusted_fits.(p), this_data);
            plot(time, model_wtpercent, 'Color', colours{line_id});
            
            line_id = line_id + 1;
            
        end
    end

    ax = gca();
    ax.FontName = font;
    ax.FontSize = 11;
    ax.XLim = [-4 22];
    ax.XTick = [-4 0:5:20];
    ax.XTickLabel{1} = 'BR';
    ax.TickLength = [0.02 0.05];
%     ax.TickDir = 'out';
    title("(" + char('a'-1+subplot_id) + ") " + p + " with adjusted boundary layer delta");

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";


%% Functions

function [time, model_wtpercent, model_conc, associated, free_plasticiser, coords] = predict(fit, data, time)
    R = 8.31446; % universal gas constant
    size = unique(data.Size);
    wc = unique(data.WaterCondition);
    temp = unique(data.Temperature);
    assert(numel(size) == 1);
    assert(numel(wc) == 1);
    assert(numel(temp) == 1);

    r0 = data.Properties.UserData.radius(size);
    if wc == "Still"
        delta = fit.delta_still;
    elseif wc == "Agitated"
        delta = fit.delta_agitated;
    else
        error("Unknown water condition");
    end

    temp_correction = exp(fit.D_plastic_EA / (R * (25.5 + 273.15))) * exp(-fit.D_plastic_EA / (R * (temp + 273.15)));
    D_plastic = fit.D_plastic_25C * temp_correction;
    p = fit.plasticiser;

    associated = data.Properties.UserData.associated.(p).("temp_"+floor(temp)).(string(wc));
    p0 = data.InternalPlasticiserConc(data.Time == 0) - associated;
    assert(numel(p0) == 1);

    mdl = make_diffusion_model(200, r0, associated, D_plastic, delta);
    if nargin < 3
        time = linspace(0, 22, 300);
    end
    [model_conc, free_plasticiser] = run_model(mdl, p0, time);
    model_wtpercent = model_conc ./ (1 + model_conc) .* 100;
    coords = mdl.r;
end
