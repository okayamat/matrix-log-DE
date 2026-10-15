function [M, Sn, En] = ImprovedAlg3(A, m, max_points, erel)
    I = eye(m);

    lambdaMax = norm(A, 2);
    lambdaMin = norm(A\I, 2);
    L = max(log(lambdaMax),log(lambdaMin));
    d = acos(sqrt(2/(1 + sqrt(1 + (2*pi/L)^2))));
    q = 0.95;
    y = q*d;
    x = asinh(L / (pi*cos(y)));
    T = ((pi/2)*cos(y)/tan(y)) - L/2;
    t = T / (1 + sqrt(T));
    a = (pi/2)*sin(y) * sqrt(1 + (sinh(x) + 2*t/(pi*cos(y)))^2);
    g = cos((pi/2)*sin(y))*cos(a)*cos(y);
alpha = 2*pi*y;
    C = 2*(exp(2*q*cos(d)) + 1/g);
    nreq = max(1, ceil(log1p(C/erel) / alpha));

    K = floor((max_points - 1)/2);
    nmax = max(1, floor(K / asinh(pi*K)));
    n = min(nmax, nreq);

    N = ceil(n * asinh(pi*n));
    M = 2*N + 1;
    Sn = Sorig(A, m, n, N);
    En = C / expm1(alpha*n);

    if En > erel
        warning('ImprovedAlg3:MaxNodes', ...
            'Requested tolerance was not met: erel=%g, En=%g, M=%d.', erel, En, M);
    end
end

function Sn = Sorig(A, m, n, N)
    h  = 1 / n;
    Sn = zeros(m);
    for k=N:-1:1
        Sn = Sn + FDE_plus( A, k*h, m);
        Sn = Sn + FDE_minus(A,-k*h, m);
    end
        Sn = Sn + FDE_plus( A, 0, m);

    Sn = h * Sn;
end

function FDE = FDE_plus(A,x,m)
    I = eye(m);
    px= exp(-pi*sinh(x));
    X = (A + px*I) \ (A - I);
    weight = pi*cosh(x)*px/(1 + px);
    FDE = weight*X;
end

function FDE = FDE_minus(A,x,m)
    I = eye(m);
    px= exp(pi*sinh(x));
    X = (I + px*A) \ (A - I);
    weight = pi*cosh(x)*px/(1 + px);
    FDE = weight*X;
end
