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
data.Size = categorical(data.Size, {'<200um', '1-2mm'});
data.Water = categorical(data.Water);
data.WaterCondition = categorical(data.WaterCondition);

% sort by time 
[~, sort_key] = sort(data.Time);
data = data(sort_key, :);

% convert to g plasticiser / g plastic 
data.MeanPlasticiserConc = data.MeanPlasticiserWtPercent/100 ./ (1 - data.MeanPlasticiserWtPercent/100);

% perform fitting on <200um microplastics only 
data = data(data.Size == '<200um', :);

% read the results of the rinsing experiment 
rinsing_data = readtable('data/Rinse experiment.csv');
% convert to g plasticiser / g plastic 
rinsing_data.ConcAfterRinse = rinsing_data.ConcAfterRinseWtPercent/100 ./ (1 - rinsing_data.ConcAfterRinseWtPercent/100);

%% Define the amount of associated plasticiser (different each each temp and water condition)
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    for t = unique(data.Temperature)'
        for wc = ["Still", "Agitated"]
            this_data = data(data.Plasticiser == p & data.Temperature == t, :);
            this_data = this_data(this_data.Time >= 7, :);
            this_data = this_data(this_data.WaterCondition == wc, :);
            % DEHP agitation was not in equilibrium at 7 days, so take the
            % latest value alone
            if p == "DEHP" && wc == "Agitated"
                data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc) = min(this_data.MeanPlasticiserConc);
            else 
                data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc) = mean(this_data.MeanPlasticiserConc);
            end
        end
    end
end

%% Find the amount of surface plastic
for p = ["BPA", "BPS", "DEHT", "DEHP"]
    % find initial conc
    this_data = data(data.Plasticiser == p & data.Time == 0, :);
    start_conc = mean(this_data.MeanPlasticiserConc);

    % find rinsed conc
    this_rinsing_data = rinsing_data(rinsing_data.Plasticiser == p, :);
    rinsed_conc = this_rinsing_data.ConcAfterRinse;
    assert(numel(rinsed_conc) == 1);

    % the rinsed conc cannot be less than the total associated plasticiser 
    max_ap = 0;
    for t = unique(data.Temperature)'
        for wc = ["Still", "Agitated"]
            if ~isempty(data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc))
                max_ap = max(max_ap, data.Properties.UserData.associated.(p).("temp_"+floor(t)).(wc));
            end
        end
    end
    rinsed_conc = max(rinsed_conc, max_ap);

    % set the amount of surface plastic
    surface_plastic = max(0, start_conc - rinsed_conc);

    % save into data table
    data.Properties.UserData.surface.(p) = surface_plastic;
end

%% Set the plasticiser to curve fit
for plasticiser = ["BPA", "BPS", "DEHT", "DEHP"]

    %% Select data to use for fitting
    all_fit_data = data(data.Plasticiser == plasticiser, :);

    %% Curve fitting (25C only)

    % order of params: D_plastic_25C, D_plastic_EA, delta_still, delta_agitated
    fit_data = all_fit_data(all_fit_data.Temperature == 25.5, :);

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-2);
    lb = [1 10 10];
    ub = [5e3 1e4 1e4];
    params1 = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), 0, p(2), p(3)), numel(lb), lb, ub, options);

    [~,model_output] = evaluate_model(fit_data, plasticiser, params1(1), 0, params1(2), params1(3));
    make_plot(model_output);

    %% Curve fitting (impact of temp)

    fit_data = all_fit_data;

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-2);
    lb = [1 1 ];
    ub = [50e3 1e5];
    betas = [params1(2) params1(3)];

    params2 = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), p(2), betas(1), betas(2)), 2, lb, ub, options);

    [~,model_output] = evaluate_model(fit_data, plasticiser, params2(1), params2(2), betas(1), betas(2));
    make_plot(model_output);

    %% Curve fitting (global to fine tune the previous fit)
    fit_data = all_fit_data;

    options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100, 'UseParallel', true, 'FunctionTolerance', 1e-3);
    last_fit = [params2(1), params2(2), betas];
    lb = last_fit*0.1;
    ub = last_fit*10;
    all_params = particleswarm(@(p)evaluate_model(fit_data, plasticiser, p(1), p(2), p(3), p(4)), 4, lb, ub, options);

    [loss,model_output] = evaluate_model(data, plasticiser, all_params(1), all_params(2), all_params(3), all_params(4));
    make_plot(model_output);
    xlim([-1 25]);
    model_output.all_params = all_params;

    %%
    save(fullfile("Results", plasticiser), "model_output");

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
    markers = {'o', '^', 's', 'd'};
    
    figure();
    subplot(1,1,1);
    sgtitle(model_output.plasticiser);
    line_id = 1;

    data = model_output.data;

    % loop over each size
    for s = unique(data.Size)'
        for t = unique(data.Temperature)'
            for wc = unique(data.WaterCondition)'
                this_data = data(data.Temperature == t & data.WaterCondition == wc & data.Size == s, :);
                if ~isempty(this_data)
                    % plot 
                    errorbar(this_data.Time, this_data.MeanPlasticiserWtPercent, this_data.StdPlasticiserWtPercent, markers{line_id}, 'Color', colours{line_id}, 'MarkerEdgeColor', colours{line_id}, 'MarkerFaceColor', colours{line_id}, 'DisplayName', sprintf('Measurement: Temp = %g deg C, %s water', t, wc));
                    hold on;

                    label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));

                    plot(model_output.time, model_output.(label)*100, 'Color', colours{line_id}, 'DisplayName', sprintf('Model: Temp = %g deg C, %s water', t, wc));
                    line_id = line_id + 1;
                end
            end
        end
    end
    
    xlabel('Time (days)');
    ylabel('Plasticiser concentration (wt.%)');
    legend show;
    title(s);
    drawnow
end

