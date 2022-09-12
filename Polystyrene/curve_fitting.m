% This is the code to run the curve fitting process for the experimental
% data in the 'data' folder.
%
% System of units: g, um, days.
%
% Concentration is measured in g plasticiser / g plastic in solids, or
% g plasticiser / g water in the liquid.
%
% Diffusion coefficients are measured in um^2/day.
%

% Make sure the output directory exists 
[~,~,~] = mkdir("Results");

%% Read data
data = readtable('data/summary.csv');
% remove units from variable names
data.Properties.VariableNames = cellfun(@(s)table(split(s, '_')).Var1{1}, data.Properties.VariableNames, 'UniformOutput', false);

data.Plasticiser = categorical(data.Plasticiser);
data.Size = categorical(data.Size, {'<200um', '400um-1mm', '1-2mm'});
data.Water = categorical(data.Water);
data.WaterCondition = categorical(data.WaterCondition);

% sort by time 
[~, sort_key] = sort(data.Time);
data = data(sort_key, :);

% convert to g plasticiser / g plastic 
data.MeanPlasticiserConc = data.MeanPlasticiserWtPercent/100 ./ (1 - data.MeanPlasticiserWtPercent/100);

% specify representative radii for each size (in um)
data.Properties.UserData.radius = [136 593 1447] / 2; % convert diameter to radius, using mean values from image analysis

% read the results of the rinsing experiment
rinsing_data = readtable('data/Rinse experiment.csv');
% convert to g plasticiser / g plastic 
rinsing_data.ConcAfterRinse = rinsing_data.ConcAfterRinseWtPercent/100 ./ (1 - rinsing_data.ConcAfterRinseWtPercent/100);

%% Define the amount of associated plasticiser (different each each temp and water condition)
% Use the small microplastics only (since the larger pieces appear not to
% have reached equilibrium by the end of the experiment)
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    for t = unique(data.Temperature)'
        for wc = ["Still", "Agitated"]
            this_data = data(data.Plasticiser == p & data.Temperature == t, :);
            this_data = this_data(this_data.Time >= 7, :);
            this_data = this_data(this_data.WaterCondition == wc, :);
            this_data = this_data(this_data.Size == "<200um", :);
            if isempty(this_data)
                continue;
            end
            % DEHP was not in equilibrium at 7 days, so take the latest
            % value alone
            if p == "DEHP" || p == "DEHT"
                this_data = this_data(this_data.Time >= 21, :);
                data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc) = mean(this_data.MeanPlasticiserConc);
            else 
                data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc) = mean(this_data.MeanPlasticiserConc);
            end
        end
    end
end

%% Find the amount of surface plasticiser
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    % find initial conc
    this_data = data(data.Plasticiser == p & data.Time == 0 & data.Size == "<200um", :);
    start_conc = mean(this_data.MeanPlasticiserConc);

    % find rinsed conc
    this_rinsing_data = rinsing_data(rinsing_data.Plasticiser == p, :);
    rinsed_conc = this_rinsing_data.ConcAfterRinse;
    assert(numel(rinsed_conc) == 1);

    % the rinsed conc cannot be less than the total associated plasticiser 
    max_ap = 0;
    for t = unique(data.Temperature)'
        for wc = ["Still", "Agitated"]
            if isfield(data.Properties.UserData.associated.(p).("temp_"+floor(t)), wc)
                max_ap = max(max_ap, data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc));
            end
        end
    end
    rinsed_conc = max(rinsed_conc, max_ap);

    % set the amount of surface plasticiser
    surface_plasticiser = max(0, start_conc - rinsed_conc);

    % Scale to the other sizes. This works by assuming a constant surface
    % density (weight of plasticiser per surface area of sphere, g/um^2).
    % Apply the same surface density to a larger sphere, and convert the
    % resulting amount of surface plasticiser to an equivalent
    % concentration over the total mass of plastic. The resulting algebra
    % gives this relationship for scaling the surface plasticiser:
    r1 = data.Properties.UserData.radius(1);
    surface_plasticiser = surface_plasticiser .* r1 ./ data.Properties.UserData.radius;
    
    % save into data table
    data.Properties.UserData.surface.(p) = surface_plasticiser;
end

%% Subtract off the surface plasticiser at t=0
data.InternalPlasticiserConc = data.MeanPlasticiserConc;
data.InternalPlasticiserWtPercent = data.MeanPlasticiserWtPercent;

for p = ["BPA", "BPS", "DEHT", "DEHP"]
    for s = unique(data.Size)'
        for wc = unique(data.WaterCondition)'
            for temp = unique(data.Temperature)'
                mask = data.Plasticiser == p & data.Size == s & data.WaterCondition == wc & data.Temperature == temp & data.Time == 0;
                if ~any(mask)
                    continue;
                end
                assert(sum(mask) == 1);

                s_conc = data.Properties.UserData.surface.(p)(s);
                initial_c = data.MeanPlasticiserConc(mask) - s_conc;
                initial_c = max(initial_c, data.Properties.UserData.associated.(p).("temp_"+floor(temp)).(string(wc)));

                data.InternalPlasticiserConc(mask) = initial_c;
                c = data.InternalPlasticiserConc(mask);
                data.InternalPlasticiserWtPercent(mask) = c ./ (c + 1) * 100; % convert to wt.%
            end
        end
    end
end

%% Save the data table
save("Results\data.mat", "data");

%% Fit on smaller particles only
all_fit_data = data;
all_fit_data = all_fit_data(~(all_fit_data.Size == "1-2mm"), :);

%% Set the plasticiser to curve fit
for plasticiser = ["BPA", "BPS", "DEHT", "DEHP"]

    %% Select data to use for fitting
    plasticiser_data = all_fit_data(all_fit_data.Plasticiser == plasticiser, :);

    %% Curve fitting (25C only)

    % order of params: D_plastic_25C, D_plastic_EA, delta_still,
    % delta_agitated_factor
    fit_data = plasticiser_data(plasticiser_data.Temperature == 25.5, :);

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-2);
    lb = [1 1 1];
    ub = [1e3 1e4 1e2];
    params1 = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), 0, p(2), p(3)), numel(lb), lb, ub, options);
    [~,model_output] = evaluate_model(fit_data, plasticiser, params1(1), 0, params1(2), params1(3));
    make_plot(model_output);

    %% Curve fitting (impact of temp)

    fit_data = plasticiser_data;

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-2);
    lb = [1 1 ];
    ub = [50e3 1e4];
    betas = [params1(2) params1(3)];

    params2 = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), p(2), betas(1), betas(2)), 2, lb, ub, options);

    [~,model_output] = evaluate_model(fit_data, plasticiser, params2(1), params2(2), betas(1), betas(2));
    make_plot(model_output);

    %% Curve fitting (global to fine tune the previous fit)
    fit_data = plasticiser_data;

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-3);
    last_fit = [params2(1), params2(2), betas];
    lb = last_fit*0.1;
    lb(4) = 1;
    ub = last_fit*10;
    all_params = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), p(2), p(3), p(4)), 4, lb, ub, options);

    [loss,model_output] = evaluate_model(fit_data, plasticiser, all_params(1), all_params(2), all_params(3), all_params(4));
    make_plot(model_output);
    

    fit_result = struct();
    fit_result.plasticiser = plasticiser;
    fit_result.D_plastic_25C = all_params(1);
    fit_result.D_plastic_EA = all_params(2);
    fit_result.delta_still = all_params(3);
    fit_result.delta_agitated = all_params(3) / all_params(4);

    %%
    save(fullfile("Results", plasticiser), '-struct', "fit_result");

    %%
end

%% Functions
function make_plot(model_output)
    
    colours = {...
        [0.890625,0.1015625,0.109375], ...
        [0.21484375,0.4921875,0.71875], ...
        [0.30078125,0.68359375,0.2890625], ...
        [0.59375,0.3046875,0.63671875], ...
        [0.99609375,0.49609375,0], ...
        [0.6484375,0.3359375,0.15625], ...
        [0.96484375,0.50390625,0.74609375], ...
        [0.5,0.5,0.5], ...
        [0.8965,0.8965,0.1793]};
    markers = {'o', '^', 's', 'd', '*', 'x', '>', '<', 'p', 'v'};
    
    figure();
%     tiledlayout('flow');
    sgtitle(model_output.plasticiser);
    line_id = 1;

    data = model_output.data;

    % loop over each size
    for s = unique(data.Size)'
%         nexttile();
        for t = unique(data.Temperature)'
            for wc = unique(data.WaterCondition)'
                this_data = data(data.Temperature == t & data.WaterCondition == wc & data.Size == s, :);
                if ~isempty(this_data)
                    % plot 
                    errorbar(this_data.Time, ...
                        this_data.InternalPlasticiserWtPercent, ...
                        this_data.StdPlasticiserWtPercent, ...
                        markers{line_id}, ...
                        'Color', colours{line_id}, ...
                        'MarkerEdgeColor', colours{line_id}, ...
                        'MarkerFaceColor', colours{line_id}, ...
                        'DisplayName', sprintf('Measurement: Size=%s, Temp=%g deg C, %s water', s, t, wc));
                    hold on;

                    label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));

                    plot(model_output.time, model_output.(label)*100, 'Color', colours{line_id}, 'DisplayName', sprintf('Model: Size=%s, Temp = %g deg C, %s water', s, t, wc));
                    line_id = line_id + 1;
                end
            end
        end
%         title(s);
        xlabel('Time (days)');
        ylabel('Plasticiser concentration after rinse (wt.%)');
        legend show;
    end
    
    drawnow
end

function [loss, model_output] = evaluate_model(data, plasticiser, D_plastic_25C, D_plastic_EA, delta_still, delta_agitated_factor)
    % initialise with zero loss
    loss = 0;

    R = 8.31446; % universal gas constant
    delta_agitated = delta_still /  delta_agitated_factor;

    % filter to only the selected plasticiser
    data = data(data.Plasticiser == plasticiser, :);

    % save the model_output struct if requested
    if nargout >= 2
        model_output = struct();
        model_output.time = linspace(0, 30, 1000);
        model_output.plasticiser = plasticiser;
        model_output.data = data;
    end
    
    % loop over each size
    sizes = unique(data.Size);
    num_sizes = numel(sizes);
    for s_id = 1:num_sizes
        s = sizes(s_id);
        r0 = data.Properties.UserData.radius(s);

        temps = unique(data.Temperature);
        for t_id = 1:numel(temps)
            t = temps(t_id);
            % calculate the adjustment factor for the Arrhenius temperature
            % dependence
            temp_correction = exp(D_plastic_EA / (R * (25.5 + 273.15))) * exp(-D_plastic_EA / (R * (t + 273.15)));
            D_plastic = D_plastic_25C * temp_correction;
            for wc = unique(data.WaterCondition)'
                if wc == "Still"
                    delta = delta_still;
                elseif wc == "Agitated"
                    delta = delta_agitated;
                else
                    error("Unknown water condition");
                end
                delta = delta * temp_correction; % since delta* is proportional to D_plastic as well
                
                this_data = data(data.Temperature == t & data.WaterCondition == wc & data.Size == s, :);
                if ~isempty(this_data)
                    % evaluate the model 
                    if nargout >= 2
                        % higher quality for final plots
                        N_terms = 200;
                    else
                        % faster for initial curve fitting
                        N_terms = 100;
                    end

                    associated_plasticiser = data.Properties.UserData.associated.(plasticiser).("temp_"+floor(t)).(string(wc));
                    mdl = make_diffusion_model(N_terms, r0, associated_plasticiser, D_plastic, delta);
                    p0 = this_data.InternalPlasticiserConc(this_data.Time == 0) - associated_plasticiser;
                    if numel(p0) ~= 1
                        disp('error');
                    end
                    assert(numel(p0) == 1);
                    p = run_model(mdl, p0, this_data.Time);
                    err = abs(p' - this_data.InternalPlasticiserConc);
                    loss = loss + sum(err.^2) / numel(this_data.InternalPlasticiserConc);
                    
                    if nargout >= 2
                        p = run_model(mdl, p0, model_output.time);
                        wt_percent = p ./ (1 + p);
                        label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));
                        model_output.(label) = wt_percent;
                    end
                end
            end
        end
    end
end
