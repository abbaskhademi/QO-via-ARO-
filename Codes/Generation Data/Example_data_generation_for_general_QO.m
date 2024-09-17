% Define parameters

n = 50;  % number of variables
m= round(n/5);  % number of constraints
p_values = [0.1, 0.5, 0.9];  % percentages of positive eigenvalues
num_instances_per_group = 5;
density = 0.8;  % density for sparse matrix generation

% Loop through different p values
for group = 1:length(p_values)
    p = p_values(group);
    
    % Generate instances for each group
    for instance = 1:num_instances_per_group
        % Set seed based on group and instance number for reproducibility
        rng(22+instance);

        % Generate eigenvalues
        num_positive = round(p * n);
        eigenvals = [5 * rand(num_positive, 1); -5 * rand(n - num_positive, 1)];
        eigenvals = eigenvals(randperm(n)); % Shuffle the eigenvalues

        % Generate Q matrix
        D = diag(eigenvals);
        [R, ~] = qr(randn(n));
        Q = R * D * R';

        % Generate c vector
        c = -5 + 10 * rand(n, 1);

        % Generate A matrix
        A = full(sprand(m-1, n, density));
        ind_A = find(A);
        A(ind_A) = -5 + 10 * A(ind_A);
        A = [A; ones(1, n)];

        % Generate b vector
        b = 0.5 * A * ones(n, 1);

        % Save the instance with the new filename format
        filename = sprintf('# Problem_QO_mx%dnx%d(%.2f)_%d.mat',  m, n, p, instance);
        save(filename, 'Q', 'A', 'b', 'c');
        
        fprintf('Saved instance: %s\n', filename);
    end
end
