%% Read data
data = readtable(fullfile("Data", "data.csv"));
rinsing = readtable(fullfile("Data", "rinsing.csv"));

data.Plasticiser = categorical(data.Plasticiser);
data.Time = hours(data.Time_hours);

[~,idx] = sort(data.Time);
data = data(idx, :);

% When is steady state reached?
steady_state = struct();
steady_state.BPA = hours(300);
steady_state.DEHP = hours(12); % hours

% convert to g plasticiser / g plastic 
data.MeanPlasticiserConc = data.MeanPlasticiserWtPercent/100 ./ (1 - data.MeanPlasticiserWtPercent/100);
rinsing.MeanPlasticiserConc = rinsing.MeanPlasticiserWtPercent/100 ./ (1 - rinsing.MeanPlasticiserWtPercent/100);

% Store the t=0 (after rinsing) values
for p = unique(data.Plasticiser)'
    data.Properties.UserData.after_rinsing.(string(p)) = rinsing.MeanPlasticiserConc(rinsing.Plasticiser == p & rinsing.Condition == "After rinsing");
    data.Properties.UserData.after_rinsing_wtpc.(string(p)) = rinsing.MeanPlasticiserWtPercent(rinsing.Plasticiser == p & rinsing.Condition == "After rinsing");
    data.Properties.UserData.after_rinsing_std.(string(p)) = rinsing.StdPlasticiserWtPercent(rinsing.Plasticiser == p & rinsing.Condition == "After rinsing");
end

% Find the associated plasticiser 
for p = unique(data.Plasticiser)'
    for t = unique(data.Temperature_degC)'
        this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);
        associated = this_data.MeanPlasticiserConc(this_data.Time >= steady_state.(string(p)));
        associated = mean(associated);
        data.Properties.UserData.associated.(string(p)).("temp_"+floor(t)) = associated;
        data.Properties.UserData.associated_wtpc.(string(p)).("temp_"+floor(t)) = 100*associated ./ (1 - associated);
    end
end

% Prepare the folder for the figure outputs
[~,~,~] = mkdir("Figures");

%% Curve fit the DEHP data to an exponential

colours = [
    0.00,0.45,0.74;
    0.93,0.69,0.13;
    0.85,0.33,0.10
];

for p = ["DEHP"]
    after_rinsing_mean = data.Properties.UserData.after_rinsing_wtpc.(p);
    after_rinsing_std = data.Properties.UserData.after_rinsing_std.(p);
    
    this_data = data(data.Plasticiser == p, :);
    times = this_data.Time;
    
    associated = this_data.MeanPlasticiserWtPercent(this_data.Time >= steady_state.(p));
    associated = mean(associated);
    
    % fit the model 
    mdl = @(x,xdata) (after_rinsing_mean-associated) * exp(-xdata / x(1)) + associated;
    x = lsqcurvefit(mdl, 1, hours(this_data.Time), this_data.MeanPlasticiserWtPercent);
    fprintf("** %s time constant = %g seconds\n", p, x*60*60); % convert from hours

    % plot (minutes)
    figh = figure('Position', [300 300 453 344]);
    for t = unique(data.Temperature_degC)'
        this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);

        time = [0; minutes(this_data.Time)];
        conc = [after_rinsing_mean; this_data.MeanPlasticiserWtPercent];
        std = [after_rinsing_std; this_data.StdPlasticiserWtPercent];
        errorbar(time, conc, std, "o--", "DisplayName", sprintf("%g °C", t));
        hold on;
    end
    time = linspace(0, max(this_data.Time), 50000);
    plot(minutes(time), mdl(x, hours(time)), "k-", "DisplayName", "Model fit");

    legend("show", "Location", "SouthEast");
    ylim([8.5 17]);
    xlim([0 60])
    title(p);
    xlabel("Leaching time (minutes)");
    ylabel("Concentration of plasticiser (wt%)");

    % also make a zoomed in version
    ax = axes('Position', [0.2589 0.6541 0.6250 0.2151]);
    for t = unique(data.Temperature_degC)'
        this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);

        time = [0; days(this_data.Time)];
        conc = [after_rinsing_mean; this_data.MeanPlasticiserWtPercent];
        std = [after_rinsing_std; this_data.StdPlasticiserWtPercent];
        errorbar(time, conc, std, "o--", "DisplayName", sprintf("%g °C", t));
        hold on;
    end
    time = linspace(0, max(this_data.Time), 50000);
    plot(days(time), mdl(x, hours(time)), "k-", "DisplayName", "Model fit");
    xlabel("Leaching time (days)");
    ylabel("Concentration (wt%)");
    exportgraphics(figh, fullfile("Figures", "Model fit " + p + ".pdf"));
    exportgraphics(figh, fullfile("Figures", "Model fit " + p + ".png"), "Resolution", 600);
end


%% Fit BPA to the diffusion model 

colours = [
    0.00,0.45,0.74;
    0.93,0.69,0.13;
    0.85,0.33,0.10
];

for plasticiser = ["BPA"]
    load(fullfile("Microscope Images", "size_contribution_"+plasticiser+".mat"));

    % exclude zeros
    radii = radii(weighting > 0);
    weighting = weighting(weighting > 0);
    
    % Select the data to use for fitting
    fit_data = data(data.Plasticiser == plasticiser, :);
    fit_data.Properties.UserData.radii = radii;
    fit_data.Properties.UserData.weighting = weighting;

    % Set optimisation options
    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 40, 'UseParallel', true, 'FunctionTolerance', 1e-2);
    lb = [1 1 0.1];
    ub = [50e3 1e6 1e2];
    params1 = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), p(2), p(3)), numel(lb), lb, ub, options);

    % Run the final model
    [~,model_output] = evaluate_model(fit_data, plasticiser, params1(1), params1(2), params1(3));

    % Make the plot
    figh = figure('Position', [300 300 453 344]);
    colour_id = 1;
    for t = unique(data.Temperature_degC)'
        this_data = data(data.Plasticiser == plasticiser & data.Temperature_degC == t, :);

        % plot the data 
        time = [0; days(this_data.Time)];
        conc = [this_data.Properties.UserData.after_rinsing_wtpc.(plasticiser); this_data.MeanPlasticiserWtPercent];
        std = [this_data.Properties.UserData.after_rinsing_std.(plasticiser); this_data.StdPlasticiserWtPercent];
        errorbar(time, conc, std, "o--", "Color", colours(colour_id, :), "DisplayName", sprintf("%g °C", t));
        hold on;
            
        % plot the model 
        plot(model_output.time, model_output.("Temp"+t), "Color", colours(colour_id, :), "DisplayName", "Model fit");
        colour_id = colour_id + 1;
    end
    xlim tight;
    legend("show", "Location", "Best");
    title(plasticiser);
    xlabel("Leaching time (days)");
    ylabel("Concentration of plasticiser (wt%)");
    exportgraphics(figh, fullfile("Figures", "Model fit " + plasticiser + ".pdf"));
    exportgraphics(figh, fullfile("Figures", "Model fit " + plasticiser + ".png"), "Resolution", 600);

    % Plot the leaching of each size
    figh = figure('Position', [1000 899 691 439]);
    to_plot = model_output.Temp26_raw;
    layout = tiledlayout(3,5);
    
    for i = 1:size(to_plot, 1)
        nexttile();
        h = plot(model_output.time, 1000*to_plot(i, :));
        title(sprintf("r = %.0f um", radii(i)));
        ax = ancestor(h, 'axes');
        ax.YAxis.Exponent = 0;
%         ytickformat('%.0f');
    end
    xlabel(layout, "Leaching time (days)");
    ylabel(layout, "Free plasticiser (mg plasticiser / g plastic)");
%     linkaxes(layout.Children, 'xy');
    sgtitle(plasticiser + " at 26 °C");
    exportgraphics(figh, fullfile("Figures", "Size resolved leaching " + plasticiser + ".pdf"));
    exportgraphics(figh, fullfile("Figures", "Size resolved leaching " + plasticiser + ".png"), "Resolution", 600);

    % Print the curve fit results 
    fprintf("** Curve fit results for %s\n", plasticiser);
    fprintf("  Diffusion coefficient at 25.5C: %.4g cm^2/s\n", params1(1) * (10^-4)^2 / 86400); % convert from um^2/days 
    fprintf("  Diffusion coefficient activation energy: %g kJ/mol\n", params1(2) / 1000); % convert from J/mol
    fprintf("  Boundary layer coefficient delta* = %g mm\n", params1(3) * 1e-6 * 1e3); % convert from um to mm


end


%% Local functions

function [loss, model_output] = evaluate_model(data, plasticiser, D_plastic_25C, D_plastic_EA, delta)
    loss = 0; % initialise with zero loss
    R = 8.31446; % universal gas constant

    % save the model_output struct if requested
    if nargout >= 2
        model_output = struct();
        model_output.time = linspace(0, 21, 1000);
        model_output.plasticiser = plasticiser;
        model_output.data = data;
    end

    % get the distribution of particle radii
    r0 = data.Properties.UserData.radii;
    r0_weights = data.Properties.UserData.weighting;
    
    % loop over each temp
    temps = unique(data.Temperature_degC);
    for t_id = 1:numel(temps)
        t = temps(t_id);
        % calculate the adjustment factor for the Arrhenius temperature
        % dependence
        temp_correction = exp(D_plastic_EA / (R * (25.5 + 273.15))) * exp(-D_plastic_EA / (R * (t + 273.15)));
        D_plastic = D_plastic_25C * temp_correction;

        this_data = data(data.Temperature_degC == t, :);
        if ~isempty(this_data)
            % evaluate the model 
            if nargout >= 2
                % higher quality for final plots
                N_terms = 200;
            else
                % faster for initial curve fitting
                N_terms = 100;
            end

            % what is the initial condition?
            associated_plasticiser = data.Properties.UserData.associated.(plasticiser).("temp_"+floor(t));
            p0 = data.Properties.UserData.after_rinsing.(plasticiser) - associated_plasticiser;

            % create a model for each radius and accumulate the result
            p = zeros(1, numel(this_data.Time));
            for r0_idx = 1:numel(r0)
                mdl = make_diffusion_model(N_terms, r0(r0_idx), 0, 0, D_plastic, delta);                
                p = p + r0_weights(r0_idx) .* run_model(mdl, p0, days(this_data.Time));
            end
            p = p + associated_plasticiser;

            % compute error in wtpc
            p_wtpc = 100 * p ./ (1 + p);
            
            err = abs(p_wtpc' - this_data.MeanPlasticiserWtPercent);
            loss = loss + sum(err.^2) / numel(this_data.MeanPlasticiserConc);
            
            % generate output for plotting
            if nargout >= 2
                label = genvarname("Temp" + t);
                p = zeros(1, numel(model_output.time));
                model_output.(label+"_raw") = zeros(numel(r0), numel(model_output.time));
                for r0_idx = 1:numel(r0)
                    mdl = make_diffusion_model(N_terms, r0(r0_idx), 0, 0, D_plastic, delta);    
                    model_output.(label+"_raw")(r0_idx, :) = r0_weights(r0_idx) .* run_model(mdl, p0, model_output.time);
                    p = p + model_output.(label+"_raw")(r0_idx, :) ;
                end
                p = p + associated_plasticiser;
                wt_percent = 100 * p ./ (1 + p);
                model_output.(label) = wt_percent;
            end
        end
    end
end
