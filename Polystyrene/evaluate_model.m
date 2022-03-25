function [loss, model_output] = evaluate_model(data, plasticiser, D_plastic_25C, D_plastic_EA, delta_still, delta_agitated)
    % initialise with zero loss
    loss = 0;

    R = 8.31446; % universal gas constant

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
        if s == "<200um"
            r0 = 100; % radius in um 
            loss_weight = 3;
        elseif s == "1-2mm"
            r0 = 1500/2; 
            loss_weight = 1;
        else
            error("Unknown size");
        end

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

                    surface_plasticiser = data.Properties.UserData.surface.(plasticiser);
                    associated_plasticiser = data.Properties.UserData.associated.(plasticiser).("temp_"+floor(t)).(string(wc));
                    mdl = make_diffusion_model(N_terms, r0, associated_plasticiser, surface_plasticiser, D_plastic, delta);
                    p0 = this_data.MeanPlasticiserConc(this_data.Time == 0) - associated_plasticiser - surface_plasticiser;
                    p = run_model(mdl, p0, this_data.Time);
                    p(this_data.Time == 0) = p(this_data.Time == 0);
                    wt_percent = p ./ (1 + p);
                    err = abs((wt_percent*100)' - this_data.MeanPlasticiserWtPercent);
                    loss = loss + loss_weight*sum(err.^2) / numel(this_data.MeanPlasticiserWtPercent);
                    
                    if nargout >= 2
                        p = run_model(mdl, p0, model_output.time);
                        p(model_output.time == 0) = p(model_output.time == 0);
                        wt_percent = p ./ (1 + p);
                        label = genvarname("Size" + string(s) + "_Temp" + t + "_" + string(wc));
                        model_output.(label) = wt_percent;
                    end
                end
            end
        end
    end
end
