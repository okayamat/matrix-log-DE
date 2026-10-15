clear all

max_points = 5*2.^(0:7);

for m = 50*2.^(0:5)
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

%%%%%%%%%% Algorithm 1 %%%%%%%%%%
    eps = 2^(-53);
    tau = 10^(-11);
    n0  = 2;

    time_history = zeros(numel(max_points), 1);
    err_history = zeros(numel(max_points), 1);
    for i = 1:numel(max_points)
      nmax = max_points(i);
      warning_state = warning('query', 'OriginalAlg1:MaxNodes');
      warning('off', 'OriginalAlg1:MaxNodes');
      warning_cleanup = onCleanup(@() warning(warning_state.state, warning_state.identifier));
      t = timeit(@() OriginalAlg1(A, m, n0, nmax, eps, tau));
      time_history(i) = t;
      [M_history, S_history] = OriginalAlg1(A, m, n0, nmax, eps, tau);
      clear warning_cleanup;
      err_history(i) = norm(logA - S_history{end}, 2) / logA2norm;
    end

    data = [time_history(:), err_history(:)];
    filename = sprintf('%d/Alg1-TimeError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Algorithm 2 %%%%%%%%%%
    eabs = 10^(-13);
    erel = 10^(-11);

    time_history = zeros(numel(max_points), 1);
    err_history = zeros(numel(max_points), 1);
    for i = 1:numel(max_points)
      nmax = max_points(i);
      warning_state = warning('query', 'ImprovedAlg2:MaxNodes');
      warning('off', 'ImprovedAlg2:MaxNodes');
      warning_cleanup = onCleanup(@() warning(warning_state.state, warning_state.identifier));
      t = timeit(@() ImprovedAlg2(A, m, nmax, eabs, erel));
      time_history(i) = t;
      [M_history, S_history] = ImprovedAlg2(A, m, nmax, eabs, erel);
      clear warning_cleanup;
      err_history(i) = norm(logA - S_history{end}, 2) / logA2norm;
    end

    data = [time_history(:), err_history(:)];
    filename = sprintf('%d/Alg2-TimeError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Algorithm 3 %%%%%%%%%%
    time_history = zeros(numel(max_points), 1);
    err_history = zeros(numel(max_points), 1);
    for i = 1:numel(max_points)
      nmax = max_points(i);
      warning_state = warning('query', 'ImprovedAlg3:MaxNodes');
      warning('off', 'ImprovedAlg3:MaxNodes');
      warning_cleanup = onCleanup(@() warning(warning_state.state, warning_state.identifier));
      t = timeit(@() ImprovedAlg3(A, m, nmax, erel));
      time_history(i) = t;
      [M, Sn, En] = ImprovedAlg3(A, m, nmax, erel);
      clear warning_cleanup;
      err_history(i) = norm(logA - Sn, 2) / logA2norm;
    end

    data = [time_history(:), err_history(:)];
    filename = sprintf('%d/Alg3-TimeError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
end
