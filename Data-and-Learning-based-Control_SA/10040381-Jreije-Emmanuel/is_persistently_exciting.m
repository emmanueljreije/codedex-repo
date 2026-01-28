function [flag, rH, nRows] = is_persistently_exciting(u, Lorder)
% u: m × N input data
% Lorder: order to check (L+n)

    % your Hankel builder expects z = N × m
    z = u.';    % transpose: now N × m

    % build block Hankel matrix using function H
    Hmat = H(Lorder, z);

    % rank test
    rH    = rank(Hmat);
    nRows = size(Hmat,1);

    flag = (rH == nRows);
end