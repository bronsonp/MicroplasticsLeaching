% This code makes the plots that are used in the paper. First run the
% `curve_fitting` script to produce the results, then run this script to
% produce the figures.

%% Load data 
results = struct();
for p = ["BPS", "BPA", "DEHT", "DEHP"]
    tmp = load(fullfile("Results", p+".mat"));
    results.(p) = tmp.model_output;
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

%% Figure 2 - impact of temperature

figh=figure();
s = "<200um";
wc = "Still";
subplot_id = 1;
tiles=tiledlayout(2,2,'TileSpacing', 'compact');
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    data = results.(p).data;

    nexttile();
    line_id = 1;
    lg = [];
    
    for t = unique(data.Temperature)'
        this_data = data(data.Temperature == t & data.WaterCondition == wc & data.Size == s, :);
        if ~isempty(this_data)
            % plot
            lg(end+1) = errorbar(this_data.Time, this_data.MeanPlasticiserWtPercent, this_data.StdPlasticiserWtPercent, markers{line_id}, 'Color', colours{line_id}, 'MarkerEdgeColor', colours{line_id}, 'MarkerFaceColor', colours{line_id}, 'DisplayName', sprintf('%g °C', t));
            hold on;
        
            label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));
        
            plot(results.(p).time, results.(p).(label)*100, 'Color', colours{line_id}, 'HandleVisibility','off');
            line_id = line_id + 1;
            p0 = this_data.MeanPlasticiserConc(this_data.Time==0);
        end
    end

    sp = data.Properties.UserData.surface.(p);
    after_rise = p0 - sp;
    after_rise = after_rise / (after_rise + 1) * 100; % convert to wt.%
    lg(end+1) = plot(get(gca(), 'xlim'), [1 1]*after_rise, 'k--', 'DisplayName', 'After 60s rinse');

    set(gca, 'FontName', font);
    
    xlim([-1 22]);
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";

%% Figure 3 - impact of agitation
figure();
s = "<200um";
t = 25.5;
subplot_id = 1;
tiles=tiledlayout(2,2,'TileSpacing', 'compact');
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    data = results.(p).data;

    nexttile();
    lg = [];
    
    line_id = 1;
    
    for wc = ["Still", "Agitated"]
        this_data = data(data.Temperature == t & data.WaterCondition == wc & data.Size == s, :);
        if ~isempty(this_data)
            % plot
            lg(end+1) = errorbar(this_data.Time, this_data.MeanPlasticiserWtPercent, this_data.StdPlasticiserWtPercent, markers{line_id*2}, 'Color', colours{line_id*2}, 'MarkerEdgeColor', colours{line_id*2}, 'MarkerFaceColor', colours{line_id*2}, 'DisplayName', sprintf('%s', wc));
            hold on;
        
            label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));
        
            plot(results.(p).time, results.(p).(label)*100, 'Color', colours{line_id*2}, 'HandleVisibility','off');
            line_id = line_id + 1;
            p0 = this_data.MeanPlasticiserConc(this_data.Time==0);
        end
    end

    set(gca, 'FontName', font);
    
    sp = data.Properties.UserData.surface.(p);
    after_rise = p0 - sp;
    after_rise = after_rise / (after_rise + 1) * 100; % convert to wt.%
    lg(end+1) = plot(get(gca(), 'xlim'), [1 1]*after_rise, 'k--', 'DisplayName', 'After 60s rinse');

    xlim([-1 22]);
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end

xlabel(tiles, 'Leaching time (days)');
ylabel(tiles, 'Plasticiser conc. (wt%)');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";

%% Table 1 

clc;
R = 8.31446; % universal gas constant
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    disp('---');
    disp(p);

    D_25C = results.(p).all_params(1);
    Ea = results.(p).all_params(2);
    delta_still = results.(p).all_params(3);
    delta_agitated = results.(p).all_params(4);

    fprintf('Diffusion coeff at 25.5 C: %.3g cm^2\n', D_25C * (10^-4)^2 / 86400);
    fprintf('Ea: %g kJ/mol\n', Ea/1000);
    fprintf('delta*_still: %.3e m\n', delta_still * 1e-6);
    fprintf('delta*_agitated: %.3e m\n', delta_agitated * 1e-6);
end

%% Evaluate the spatial behaviour of each plasticiser

R = 8.31446; % universal gas constant
s = "<200um";
temp = 25.5;
wc = "Still";

figure();
tiles=tiledlayout(2,2,'TileSpacing', 'compact');
subplot_id = 1;
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    % find the amount of time that it takes for the free plasticiser to
    % decay by half
    t = results.(p).time;
    y = results.(p).(genvarname("Size" + string(s) + "_Temp" + string(temp) + "_" + string(wc)));
    y = y ./ (1 - y); % convert to g/g
    sp = results.(p).data.Properties.UserData.surface.(p);
    y(1) = y(1) - sp; % subtract surface plasticiser
    associated_plasticiser = results.(p).data.Properties.UserData.associated.(p).("temp_"+floor(temp)).(string(wc));
    y_free = y - associated_plasticiser;
    y_target = y_free(1) * 0.01;
    mask = y_free >= y_target;
    mask(find(~mask, 1, 'first')) = true; % extend the mask by one more point
    t = t(mask);
    y_free = y_free(mask);
    t_target = interp1(y_free, t, y_target);

    % collect model parameters
    N_terms = 500;
    r0 = 100; % radius in um
    D_plastic_25C = results.(p).all_params(1);
    D_plastic_EA = results.(p).all_params(2);
    temp_correction = exp(D_plastic_EA / (R * (25.5 + 273.15))) * exp(-D_plastic_EA / (R * (temp + 273.15)));
    D_plastic = D_plastic_25C * temp_correction;
    if wc == "Still"
        delta = results.(p).all_params(3);
    else
        delta = results.(p).all_params(4);
    end

    % evaluate the model at various times
    times = t_target .* [0 0.01 0.05 0.1:0.1:1];
        mdl = make_diffusion_model(N_terms, r0, 0, 0, D_plastic, delta);
    [avg_p, free_plasticiser] = run_model(mdl, 1, times);
    
    % make the plot
    nexttile
    lg = [];
    for i = 1:numel(times)
        hsv = [(i-1)/(numel(times)) 0.9 0.8];
        rgb = hsv2rgb(hsv);
        lg(end+1) = plot(mdl.r, free_plasticiser(:, i), 'Color', rgb, 'LineWidth', 1, 'DisplayName', sprintf('t / t_{leach} = %g', times(i)/t_target));
        hold on;
    end
    title("(" + char('a'-1+subplot_id) + ") " + p);

    subplot_id = subplot_id + 1;
end
xlabel(tiles, 'Radial position (μm)');
ylabel(tiles, 'Relative concentration of plasticiser, p / p_0');
l = legend(nexttile(2), lg);
l.Location = "northeastoutside";





