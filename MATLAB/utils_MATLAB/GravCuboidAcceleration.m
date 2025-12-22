classdef GravCuboidAcceleration
% Gravitational acceleration of a homogeneous cuboid
methods(Static)

%% -----------------------------------------------
% Safe math 
%% -----------------------------------------------
function out = sp_log(a)
    if a <= 0
        out = 0;
    else
        out = log(a);
    end
end

function out = sp_atan(num, den)
    if num == 0 || den == 0
        out = 0;
    else
        out = atan(num/den);
    end
end


%% -----------------------------------------------
% Vertex coordinates & signs
%% -----------------------------------------------
function S = vertex_signs()
    S = zeros(8,1);
    idx = 1;
    for ix = 0:1
        for iy = 0:1
            for iz = 0:1
                S(idx) = (-1)^(ix+iy+iz);
                idx = idx + 1;
            end
        end
    end
end

function [X,Y,Z,S] = vertex_coordinates(x,y,z,L,B,D)
    Xb = [ L - x ; -(L + x) ];
    Yb = [ B - y ; -(B + y) ];
    Zb = [ D - z ; -(D + z) ];

    X = repelem(Xb,4);
    Y = repmat(repelem(Yb,2),2,1);
    Z = repmat(Zb,4,1);

    S = GravCuboidAcceleration.vertex_signs();
end


%% -----------------------------------------------
% Core kernel g_core(j,k,l)
%% -----------------------------------------------
function out = g_core(xj,xk,xl,Xc,Yc,Zc)
    out = zeros(8,1);
    for i = 1:8
        X = Xc(i); Y = Yc(i); Z = Zc(i);
        r = sqrt(X*X + Y*Y + Z*Z);

        term1 = xj(i) * GravCuboidAcceleration.sp_log(xk(i) + r);
        term2 = xk(i) * GravCuboidAcceleration.sp_log(xj(i) + r);

        num = xj(i)*xk(i);
        den = xl(i)*r;

        term3 = -xl(i) * GravCuboidAcceleration.sp_atan(num, den);

        out(i) = term1 + term2 + term3;
    end
end


%% -----------------------------------------------
% Acceleration components
%% -----------------------------------------------
function gx = accel_gx(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidAcceleration.vertex_coordinates(x,y,z,L,B,D);
    kernel = GravCuboidAcceleration.g_core(Y,Z,X, X,Y,Z);
    gx = -G * rho * sum(S .* kernel);
end

function gy = accel_gy(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidAcceleration.vertex_coordinates(x,y,z,L,B,D);
    kernel = GravCuboidAcceleration.g_core(X,Z,Y, X,Y,Z);
    gy = -G * rho * sum(S .* kernel);
end

function gz = accel_gz(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidAcceleration.vertex_coordinates(x,y,z,L,B,D);
    kernel = GravCuboidAcceleration.g_core(X,Y,Z, X,Y,Z);
    gz = -G * rho * sum(S .* kernel);
end


%% -----------------------------------------------
% Combined acceleration vector
%% -----------------------------------------------
function [gx,gy,gz] = accel_point(x,y,z,L,B,D,rho,G)
    gx = GravCuboidAcceleration.accel_gx(x,y,z,L,B,D,rho,G);
    gy = GravCuboidAcceleration.accel_gy(x,y,z,L,B,D,rho,G);
    gz = GravCuboidAcceleration.accel_gz(x,y,z,L,B,D,rho,G);
end


%% -----------------------------------------------
% Batch (optional parallel)
%% -----------------------------------------------
function [GX,GY,GZ] = accel_batch(points,L,B,D,rho,G,threshold)
    if nargin < 7
        threshold = 100;
    end

    pts = double(points);
    N = size(pts,1);

    GX = zeros(N,1);
    GY = zeros(N,1);
    GZ = zeros(N,1);

    if N > threshold
        parfor i = 1:N
            [GX(i),GY(i),GZ(i)] = GravCuboidAcceleration.accel_point( ...
                pts(i,1), pts(i,2), pts(i,3), L,B,D,rho,G );
        end
    else
        for i = 1:N
            [GX(i),GY(i),GZ(i)] = GravCuboidAcceleration.accel_point( ...
                pts(i,1), pts(i,2), pts(i,3), L,B,D,rho,G );
        end
    end
end

end % methods
end % classdef
