m = 100;
rng(739,'twister');
[Q,R] = qr(randn(m)); % QR factorization from randam
s = sign(diag(R));
s(s == 0) = 1;        % if s(k)==0; s(k)=1; end
Q = Q * diag(s);      % Obtain the Haar-uniform matrix
writematrix(Q, 'Qdata.csv');
