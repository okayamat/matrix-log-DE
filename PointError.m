clear all

max_points = 640;

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

    [M_history, S_history] = OriginalAlg1(A, m, n0, max_points, eps, tau);

    err_history = zeros(numel(S_history), 1);
    for k = 1:numel(S_history)
      err_history(k) = norm(logA - S_history{k}, 2) / logA2norm;
    end

    data = [M_history(:), err_history(:)];
    filename = sprintf('%d/Alg1-PointError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Algorithm 2 %%%%%%%%%%
    eabs = 10^(-13);
    erel = 10^(-11);

    [M_history, S_history] = ImprovedAlg2(A, m, max_points, eabs, erel);

    err_history = zeros(numel(S_history), 1);
    for k = 1:numel(S_history)
      err_history(k) = norm(logA - S_history{k}, 2) / logA2norm;
    end

    data = [M_history(:), err_history(:)];
    filename = sprintf('%d/Alg2-PointError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Algorithm 3 %%%%%%%%%%
    eabs = 10^(-13);
    erel = 10^(-11);

    [M, Sn, En] = ImprovedAlg3(A, m, max_points, erel);

    err_history = norm(logA - Sn, 2) / logA2norm;

    data = [M, err_history];
    filename = sprintf('%d/Alg3-PointError-1e%d.dat', m, round(log10(kappa)));
    writematrix(data, filename, 'Delimiter', 'tab');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
end
