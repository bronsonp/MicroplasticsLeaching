function [avg_p, p] = run_model(mdl, p0, t)
% Run the given model for the specified initial condition, returning
% results for the given time steps. The `mdl` structure comes from
% `make_diffusion_model`. 

    % r=0 boundary condition requires first two elements to be equal.    
    p0 = ones(mdl.N, 1) * p0;
    
    % preallocate
    p = zeros(numel(p0), numel(t));
    
    % calculate
    [V,D] = eig(full(mdl.M));
    c = V \ p0(2:end-1); 
    for i = 1:numel(t)
        p(2:end-1, i) = sum(V .* exp(diag(D)*t(i))' .* c', 2);
    end
    
    % r=0 boundary condition requires first two elements to be equal
    p(1,:) = p(2,:);
    
    % r=N-1 boundary condition is as follows
    p(end, :) = mdl.last_element_scaling * p(end-1, :);

    % Now integrate to get average concentration
    avg_p = 3 / mdl.r(end)^3 * trapz(mdl.r, reshape(mdl.r, [], 1).^2 .* (p + mdl.associated_plasticiser), 1);

    % Add in the surface plasticiser
    avg_p(t == 0) = avg_p(t == 0) + mdl.surface_plasticiser;
end
