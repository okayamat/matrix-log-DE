function [M_history, S_history] = ImprovedAlg2(A, max_points, eabs, erel)
    m = length(A);
    I = eye(m);

    n = 1;
    N = ceil(n * asinh(pi*n));
    M = 2*N + 1;
    Sn = Sorig(A, m, n, N);
    M_history(1,1) = M;
    S_history{1,1} = Sn;

    converged = false;
    NN = ceil(2*n * asinh(2*pi*n));
    while (2*NN + 1 <= max_points)
        S2n = 0.5*Sn + Sdiff(A, m, n, N, NN);
        n = 2*n;
        N = NN;
        M = 2*N + 1;
        M_history(end+1,1) = M;
        S_history{end+1,1} = S2n;
        entrywise_ok = abs(S2n - Sn) <= eabs + erel * abs(S2n);
        if all(entrywise_ok(:))
            converged = true;
            break;
        end
        Sn = S2n;
        NN = ceil(2*n * asinh(2*pi*n));
    end
    if ~converged
        warning('ImprovedAlg2:MaxNodes', ...
            'maximum node count %d; stopping criterion not met by M=%d.', max_points, M);
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

function Sdiff = Sdiff(A, m, n, N, NN)
    h  = 1 / (2*n);
    Sdiff = 0;
    for k=NN:-1:2*N+1
        Sdiff = Sdiff + FDE_plus( A, k*h, m);
        Sdiff = Sdiff + FDE_minus(A,-k*h, m);
    end
    for j=N:-1:1
        Sdiff = Sdiff + FDE_plus( A,(2*j-1)*h,m);
        Sdiff = Sdiff + FDE_minus(A,(-2*(j-1)-1)*h,m);
    end

    Sdiff = h * Sdiff;
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
