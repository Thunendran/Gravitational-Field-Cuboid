classdef GravCuboidPotential

methods(Static)

% ============================================================
% SAFE MATH
% ============================================================
function out = sp_log(argument)
    if argument <= 0
        out = 0;
    else
        out = log(argument);
    end
end

function out = sp_arctan(numer, denom)
    if numer == 0 || denom == 0
        out = 0;
    else
        out = atan(numer / denom);
    end
end

% ============================================================
% CORE HELPERS
% ============================================================
function R = R3(a,b,c)
    R = sqrt(a*a + b*b + c*c);
end

function [val, R] = log_term(a,b,c,X)
    R = GravCuboidPotential.R3(a,b,c);
    val = GravCuboidPotential.sp_log(X + R);
end

function [val, R] = atan_term(a,b,c,X)
    R = GravCuboidPotential.R3(a,b,c);
    val = GravCuboidPotential.sp_arctan(b*c, X*R);
end

% ============================================================
% FACE CONTRIBUTIONS
% ============================================================
function V = Front(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X2,X4,X6,-X2);
    [log2,~] = GravCuboidPotential.log_term(X2,X4,X5,-X2);
    [log3,~] = GravCuboidPotential.log_term(X2,X3,X5,-X2);
    [log4,~] = GravCuboidPotential.log_term(X2,X3,X6,-X2);

    [atan1,~] = GravCuboidPotential.atan_term(X2,X4,X6,X2);
    [atan2,~] = GravCuboidPotential.atan_term(X2,X4,X5,X2);
    [atan3,~] = GravCuboidPotential.atan_term(X2,X3,X5,X2);
    [atan4,~] = GravCuboidPotential.atan_term(X2,X3,X6,X2);

    F_L = -( ...
        X6*X4*log1 + ...
        X4*X5*log2 + ...
        X5*X3*log3 + ...
        X3*X6*log4 );

    F_T = -(X2*X2/2)*(atan1+atan2+atan3+atan4);

    V = F_L + F_T;
end


function V = Back(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X1,X4,X5,X1);
    [log2,~] = GravCuboidPotential.log_term(X1,X4,X6,X1);
    [log3,~] = GravCuboidPotential.log_term(X1,X3,X5,X1);
    [log4,~] = GravCuboidPotential.log_term(X1,X3,X6,X1);

    [atan1,~] = GravCuboidPotential.atan_term(X1,X4,X5,X1);
    [atan2,~] = GravCuboidPotential.atan_term(X1,X4,X6,X1);
    [atan3,~] = GravCuboidPotential.atan_term(X1,X3,X5,X1);
    [atan4,~] = GravCuboidPotential.atan_term(X1,X3,X6,X1);

    B_L =  X4*X5*log1 + X6*X4*log2 + X5*X3*log3 + X3*X6*log4;
    B_T = -(X1*X1/2)*(atan1+atan2+atan3+atan4);

    V = B_L + B_T;
end


function V = Right(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X1,X3,X6,X3);
    [log2,~] = GravCuboidPotential.log_term(X1,X3,X5,X3);
    [log3,~] = GravCuboidPotential.log_term(X2,X3,X5,X3);
    [log4,~] = GravCuboidPotential.log_term(X2,X3,X6,X3);

    [atan1,~] = GravCuboidPotential.atan_term(X3,X1,X6,X3);
    [atan2,~] = GravCuboidPotential.atan_term(X3,X1,X5,X3);
    [atan3,~] = GravCuboidPotential.atan_term(X3,X2,X5,X3);
    [atan4,~] = GravCuboidPotential.atan_term(X3,X2,X6,X3);

    R_L =  X6*X1*log1 + X1*X5*log2 + X5*X2*log3 + X2*X6*log4;
    R_T = -(X3*X3/2)*(atan1+atan2+atan3+atan4);

    V = R_L + R_T;
end


function V = Left(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X1,X4,X6,-X4);
    [log2,~] = GravCuboidPotential.log_term(X1,X4,X5,-X4);
    [log3,~] = GravCuboidPotential.log_term(X2,X4,X5,-X4);
    [log4,~] = GravCuboidPotential.log_term(X2,X4,X6,-X4);

    [atan1,~] = GravCuboidPotential.atan_term(X4,X1,X6,X4);
    [atan2,~] = GravCuboidPotential.atan_term(X4,X1,X5,X4);
    [atan3,~] = GravCuboidPotential.atan_term(X4,X2,X5,X4);
    [atan4,~] = GravCuboidPotential.atan_term(X4,X2,X6,X4);

    L_L = -( ...
        X6*X1*log1 + ...
        X1*X5*log2 + ...
        X5*X2*log3 + ...
        X2*X6*log4 );

    L_T = -(X4*X4/2)*(atan1+atan2+atan3+atan4);

    V = L_L + L_T;
end


function V = Top(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X1,X4,X6,-X6);
    [log2,~] = GravCuboidPotential.log_term(X2,X4,X6,-X6);
    [log3,~] = GravCuboidPotential.log_term(X2,X3,X6,-X6);
    [log4,~] = GravCuboidPotential.log_term(X1,X3,X6,-X6);

    [atan1,~] = GravCuboidPotential.atan_term(X6,X1,X4,X6);
    [atan2,~] = GravCuboidPotential.atan_term(X6,X4,X2,X6);
    [atan3,~] = GravCuboidPotential.atan_term(X6,X2,X3,X6);
    [atan4,~] = GravCuboidPotential.atan_term(X6,X3,X1,X6);

    T_L = -( ...
        X1*X4*log1 + ...
        X4*X2*log2 + ...
        X2*X3*log3 + ...
        X3*X1*log4 );

    T_T = -(X6*X6/2)*(atan1+atan2+atan3+atan4);

    V = T_L + T_T;
end


function V = Bottom(X1,X2,X3,X4,X5,X6)

    [log1,~] = GravCuboidPotential.log_term(X1,X4,X5,X5);
    [log2,~] = GravCuboidPotential.log_term(X2,X4,X5,X5);
    [log3,~] = GravCuboidPotential.log_term(X2,X3,X5,X5);
    [log4,~] = GravCuboidPotential.log_term(X1,X3,X5,X5);

    [atan1,~] = GravCuboidPotential.atan_term(X5,X1,X4,X5);
    [atan2,~] = GravCuboidPotential.atan_term(X5,X2,X4,X5);
    [atan3,~] = GravCuboidPotential.atan_term(X5,X2,X3,X5);
    [atan4,~] = GravCuboidPotential.atan_term(X5,X3,X1,X5);

    B_L =  X1*X4*log1 + X4*X2*log2 + X2*X3*log3 + X3*X1*log4;
    B_T = -(X5*X5/2)*(atan1+atan2+atan3+atan4);

    V = B_L + B_T;
end

% ============================================================
% SINGLE POINT POTENTIAL (exact Python translation)
% ============================================================
function V = potential_single(x,y,z,L,B,D,rho,G)

    X1 = (L - x);
    X2 = (L + x);
    X3 = (B - y);
    X4 = (B + y);
    X5 = (D - z);
    X6 = (D + z);

    faces = ...
        GravCuboidPotential.Front (X1,X2,X3,X4,X5,X6) + ...
        GravCuboidPotential.Back  (X1,X2,X3,X4,X5,X6) + ...
        GravCuboidPotential.Right (X1,X2,X3,X4,X5,X6) + ...
        GravCuboidPotential.Left  (X1,X2,X3,X4,X5,X6) + ...
        GravCuboidPotential.Top   (X1,X2,X3,X4,X5,X6) + ...
        GravCuboidPotential.Bottom(X1,X2,X3,X4,X5,X6);

    V = G * rho * faces;
end

% ============================================================
% VECTOR BATCH (with optional parallel)
% ============================================================
function V = potential_batch(points, L,B,D, rho, G, parallel_threshold)

    if nargin < 7
        parallel_threshold = 100;
    end

    pts = double(points);
    [N, cols] = size(pts);
    if cols ~= 3
        error("points must be Nx3");
    end

    V = zeros(N,1);

    if N > parallel_threshold
        fprintf("Using parallel processing...\n");
        parfor i = 1:N
            p = pts(i,:);
            V(i) = GravCuboidPotential.potential_single(p(1),p(2),p(3),L,B,D,rho,G);
        end
    else
        for i = 1:N
            p = pts(i,:);
            V(i) = GravCuboidPotential.potential_single(p(1),p(2),p(3),L,B,D,rho,G);
        end
    end
end

end % static methods

end % classdef
