function mdl = make_diffusion_model(N, r0, associated_plasticiser, D_plastic, delta)
% Initialise a diffusion model with the given parameters. The return value
% is a structure that contains precomputed matrices that can be used to
% quickly solve the diffusion equation. Pass this structure to `run_model`
% to actually use it.

    mdl = struct();
    mdl.associated_plasticiser = associated_plasticiser;
    mdl.N = N;

    % prepare the matrix
    matrixN = N - 2;
    mdl.r = linspace(0, r0, N);
    deltaR = mdl.r(2) - mdl.r(1);
    mdl.M = spalloc(matrixN, matrixN, 4 + (matrixN-2)*3);
    prefactor = D_plastic/(deltaR^2);
    
    % i=1
    mdl.M(1,1) = prefactor * (-2);
    mdl.M(1,2) = prefactor * (+2);
    
    % others
    for i = 2:(matrixN-1)
        mdl.M(i,i-1) = prefactor * (1-1/i);
        mdl.M(i,i) = prefactor * (-2);
        mdl.M(i,i+1) = prefactor * (1+1/i);
    end
    
    % i=N-2 (special case for RHS boundary)
    i = matrixN;
    mdl.M(i,i-1) = prefactor * (1-1/i);
    mdl.last_element_scaling = 1/(1 + 1/delta * deltaR);
    mdl.M(i,i) = prefactor * (-2 + (1+1/i) * mdl.last_element_scaling);
end
