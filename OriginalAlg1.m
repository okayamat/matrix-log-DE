function [M_history, S_history] = OriginalAlg1(A, m, n0, max_points, eps, tau)
    I = eye(m);

    mu    = norm(A-I,2);
    nu    = norm(A\I,2);
    rhoA  = max(abs(eig(A)));
    theta = abs(log(rhoA));
    epsmax= 3*mu*nu/(theta*(1+nu));
    epsuse= eps;
    if (epsuse >= epsmax)
        epsuse = epsmax/2;
    end

    a = min(theta*epsuse/(3*mu), 1/(2*mu));
    c = min((theta*epsuse)/(4*mu*nu), 1/(1 + 2*nu));
    %b= max([1 - (theta*epsuse)/(4*mu*nu), 2*nu/(2*nu+1)]);
    l = asinh((1/pi)*(log(a) - log1p(-a)));
    r = asinh((1/pi)*(log1p(-c) - log(c)));

    n = n0;
    Sn= Sorig(A, m, n, l, r);
    M = n+1;
    M_history(1,1) = M;
    S_history{1,1} = (A - I)*Sn;

    converged = false;
    while (2*n + 1 <= max_points)
        S2n = 0.5*Sn + Sdiff(A, m, n, l, r);
        n = 2*n;
        M = n + 1;
        M_history(end+1,1) = M;
        S_history{end+1,1} = (A - I)*S2n;
        if (norm(S2n - Sn,2)/(3*theta) <= tau)
            converged = true;
            break;
        end
        Sn = S2n;
    end
    if ~converged
        warning('OriginalAlg1:MaxNodes', ...
            'maximum node count %d; stopping criterion not met by M=%d.', max_points, M);
    end
end

function Sn = Sorig(A, m, n, l, r)
    h  = (r - l) / n;
    Sn = 0.5 * (GDE_at(A, l, m) + GDE_at(A, r, m));
    for k=1:n-1
        Sn = Sn + GDE_at(A, l+k*h, m);
    end
    Sn = h * Sn;
end

function Sdiff = Sdiff(A, m, n, l, r)
    h   = (r - l) / (2*n);
    Sdiff = 0;
    for j=1:n
        Sdiff = Sdiff + GDE_at(A, l+(2*j-1)*h, m);
    end
    Sdiff = h * Sdiff;
end

function GDE = GDE_at(A,x,m)
    I   = eye(m);
    qx  = (pi/2)*sinh(x);
    M   = (1 + tanh(qx))*(A - I) + 2*I;
    cqx = cosh(qx);
    weight = ((pi/2)*cosh(x)/cqx)/cqx;
    GDE = M \ (weight*I);
end
