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

p = "DEHP";

after_rinsing_mean = data.Properties.UserData.after_rinsing_wtpc.(p);
after_rinsing_std = data.Properties.UserData.after_rinsing_std.(p);
    
this_data = data(data.Plasticiser == p, :);
times = this_data.Time;

associated = this_data.MeanPlasticiserWtPercent(this_data.Time >= steady_state.(p));
associated = mean(associated);

% fit the model
DEHP_mdl = @(x,xdata) (after_rinsing_mean-associated) * exp(-xdata / x(1)) + associated;
DEHP_tau = lsqcurvefit(DEHP_mdl, 1, hours(this_data.Time), this_data.MeanPlasticiserWtPercent);
fprintf("** %s time constant = %g seconds\n", p, DEHP_tau*60*60); % convert from hours

%% Curve fit the BPS data to the model 

p = "BPA";
load(fullfile("Microscope Images", "size_contribution_"+p+".mat"));

% exclude zeros
radii = radii(weighting > 0);
weighting = weighting(weighting > 0);

% Select the data to use for fitting
fit_data = data(data.Plasticiser == p, :);
fit_data.Properties.UserData.radii = radii;
fit_data.Properties.UserData.weighting = weighting;

% Set optimisation options
options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 40, 'UseParallel', true, 'FunctionTolerance', 1e-2);
lb = [1 1 0.1];
ub = [50e3 1e6 1e2];
params1 = particleswarm(@(x)evaluate_model(fit_data, p, x(1), x(2), x(3)), numel(lb), lb, ub, options);

% Run the final model
[~,BPA_model_output] = evaluate_model(fit_data, p, params1(1), params1(2), params1(3));

% Print the curve fit results 
fprintf("** Curve fit results for %s\n", p);
fprintf("  Diffusion coefficient at 25.5C: %.4g cm^2/s\n", params1(1) * (10^-4)^2 / 86400); % convert from um^2/days 
fprintf("  Diffusion coefficient activation energy: %g kJ/mol\n", params1(2) / 1000); % convert from J/mol
fprintf("  Boundary layer coefficient delta* = %g mm\n", params1(3) * 1e-6 * 1e3); % convert from um to mm

%% Save 
save(fullfile("Figures", "curve_fit.mat"), "BPA_model_output", "params1", "DEHP_tau", "DEHP_mdl");

%% Load
load(fullfile("Figures", "curve_fit.mat"));

%% Figure 1

with_inset = true;

figh = figure('Position', [300 300    463   723]);
tiledlayout(2,1);

% (a) DEHP plot 
ax1 = nexttile();
p = "DEHP";
after_rinsing_mean = data.Properties.UserData.after_rinsing_wtpc.(p);
after_rinsing_std = data.Properties.UserData.after_rinsing_std.(p);
for t = unique(data.Temperature_degC)'
    this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);

    time = [0; days(this_data.Time)];
    conc = [after_rinsing_mean; this_data.MeanPlasticiserWtPercent];
    std = [after_rinsing_std; this_data.StdPlasticiserWtPercent];
    errorbar(time, conc, std, "o--", "DisplayName", sprintf("%g °C", t));
    hold on;
end
time = linspace(0, max(this_data.Time), 50000);
plot(days(time), DEHP_mdl(DEHP_tau, hours(time)), "k-", "LineWidth", 1.3, "DisplayName", "Model fit");
title("(a) " + p);
xlabel("Leaching time (days)");
ylabel("Concentration of plasticiser (wt%)");
xlim([0 7]);
ylim tight;
legend("show", "Location", "NorthEast");
ax1.TickDir = "out";

% (b) BPA plot
p = "BPA";
ax2 = nexttile();
colours = [
    0.00,0.45,0.74;
    0.93,0.69,0.13;
    0.85,0.33,0.10
];

% Make the plot
colour_id = 1;
for t = unique(data.Temperature_degC)'
    this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);

    % plot the data 
    time = [0; days(this_data.Time)];
    conc = [this_data.Properties.UserData.after_rinsing_wtpc.(p); this_data.MeanPlasticiserWtPercent];
    std = [this_data.Properties.UserData.after_rinsing_std.(p); this_data.StdPlasticiserWtPercent];
    errorbar(time, conc, std, "o--", "Color", colours(colour_id, :), "DisplayName", sprintf("%g °C", t));
    hold on;
        
    % plot the model 
    plot(BPA_model_output.time, BPA_model_output.("Temp"+t), "Color", colours(colour_id, :), "LineWidth", 1.3, "DisplayName", "Model fit");
    colour_id = colour_id + 1;
end
xlim([0 max(days(this_data.Time))]);
ylim([5.5 13.5]);
legend("show", "Location", "NorthEast");
title("(b) " + p);
xlabel("Leaching time (days)");
ylabel("Concentration of plasticiser (wt%)");
ax2.TickDir = "out";


% Make the inset
if with_inset 
    p = "DEHP";
    ax = axes('Position', [ 0.2947    0.7965    0.2303    0.1209]);
    for t = unique(data.Temperature_degC)'
        this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);
    
        time = [0; minutes(this_data.Time)];
        conc = [after_rinsing_mean; this_data.MeanPlasticiserWtPercent];
        std = [after_rinsing_std; this_data.StdPlasticiserWtPercent];
        errorbar(time, conc, std, "o--", "DisplayName", sprintf("%g °C", t));
        hold on;
    end
    time = linspace(0, max(this_data.Time), 50000);
    plot(minutes(time), DEHP_mdl(DEHP_tau, hours(time)), "k-", "LineWidth", 1.3, "DisplayName", "Model fit");
    xlabel("Time (mins)");
    ylabel("Conc (wt%)");
    xlim([0 15]);
    ylim tight;
end

if with_inset
    exportgraphics(figh, fullfile("Figures", "Figure 1 with inset.png"), "Resolution", 600);
else
    exportgraphics(figh, fullfile("Figures", "Figure 1.png"), "Resolution", 600);
end


%% DEHP 0-60 mins

figure();
p = "DEHP";
after_rinsing_mean = data.Properties.UserData.after_rinsing_wtpc.(p);
after_rinsing_std = data.Properties.UserData.after_rinsing_std.(p);
for t = unique(data.Temperature_degC)'
    this_data = data(data.Plasticiser == p & data.Temperature_degC == t, :);

    time = [0; minutes(this_data.Time)];
    conc = [after_rinsing_mean; this_data.MeanPlasticiserWtPercent];
    std = [after_rinsing_std; this_data.StdPlasticiserWtPercent];
    errorbar(time, conc, std, "o--", "DisplayName", sprintf("%g °C", t));
    hold on;
end
time = linspace(0, max(this_data.Time), 50000);
plot(minutes(time), DEHP_mdl(DEHP_tau, hours(time)), "k-", "LineWidth", 1.3, "DisplayName", "Model fit");
title("(a) " + p);
xlabel("Leaching time (days)");
ylabel("Concentration of plasticiser (wt%)");
xlim([0 60]);
ylim tight;
legend("show", "Location", "NorthEast");
ax1.TickDir = "out";

%% Figure 2

p = "BPA";
load(fullfile("Microscope Images", "size_contribution_"+p+".mat"));

% exclude zeros
radii = radii(weighting > 0);
weighting = weighting(weighting > 0);

for temp = [26 43 60]
    
    % Plot the leaching of each size
    figh = figure('Position', [1000 899 705 428]);
    layout = tiledlayout(3,5);
    
    for i = 1:size(to_plot, 1)
        nexttile();
    
        h = plot(BPA_model_output.time, 1000*BPA_model_output.("Temp"+num2str(temp)+"_raw")(i, :), "LineWidth", 1.3);
        
        title(sprintf("r = %.0f um", radii(i)));
        ax = ancestor(h, 'axes');
        ax.YAxis.Exponent = 0;
        ax.XTick = 0:20:60;
    %         ytickformat('%.0f');
    end
    xlabel(layout, "Leaching time (days)");
    ylabel(layout, "Free plasticiser (mg plasticiser / g plastic)");
    %     linkaxes(layout.Children, 'xy');
    sgtitle(p + " at " + num2str(temp) + " °C");
    exportgraphics(figh, fullfile("Figures", sprintf("Figure 2 T=%g C.png", temp)), "Resolution", 600);
end


%% Local functions

function [loss, model_output] = evaluate_model(data, plasticiser, D_plastic_25C, D_plastic_EA, delta)
    loss = 0; % initialise with zero loss
    R = 8.31446; % universal gas constant

    % save the model_output struct if requested
    if nargout >= 2
        model_output = struct();
        model_output.time = linspace(0, 60, 1000);
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
