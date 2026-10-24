clear all

max_points = 10:10:200;
m = 50;

  I = eye(m);
  dataname = sprintf('%d/Qdata.csv', m);
  Q = readmatrix(dataname);
  for kappa = [1e1 1e4 1e7]
    ratio = kappa^(1/(m-1));
    d = (kappa^(-1/2)) * ratio.^(0:m-1);
    D = diag(d);
    A = Q * D * Q';
    logA = Q * diag(log(d)) * Q';
    logA2norm = norm(logA,2); % 0.5 * log(kappa);

%%%%%%%%%% Algorithm 3 %%%%%%%%%%
    erel = 10^(-11);

    M_history = zeros(numel(max_points), 1);
    err_history = zeros(numel(max_points), 1);
    for i = 1:numel(max_points)
      nmax = max_points(i);
      warning_state = warning('query', 'ImprovedAlg3:MaxNodes');
      warning('off', 'ImprovedAlg3:MaxNodes');
      warning_cleanup = onCleanup(@() warning(warning_state.state, warning_state.identifier));
      [M, Sn, En] = ImprovedAlg3(A, nmax, erel);
      clear warning_cleanup;
      M_history(i) = M;
      err_history(i) = En;
    end

    data = [M_history(:), err_history(:)];
    filename = sprintf('%d/Alg3-ErrorBound-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
