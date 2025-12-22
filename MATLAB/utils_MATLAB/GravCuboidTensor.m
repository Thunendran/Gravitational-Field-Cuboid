classdef GravCuboidTensor
% Gravitational tensor of a homogeneous cuboid

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
% Vertex signs and coordinates
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

    S = GravCuboidTensor.vertex_signs();
end


%% -----------------------------------------------
% Core tensor helpers
%% -----------------------------------------------
function vals = tensor_ln(xrem, r)
    N = length(r);
    vals = zeros(N,1);
    for i = 1:N
        vals(i) = GravCuboidTensor.sp_log(xrem(i) + r(i));
    end
end

function vals = tensor_atan(xi,xj,xk,r)
    N = length(r);
    vals = zeros(N,1);
    for i = 1:N
        vals(i) = GravCuboidTensor.sp_atan(xj(i)*xk(i), xi(i)*r(i));
    end
end


%% -----------------------------------------------
% Diagonal components: Vxx, Vyy, Vzz
%% -----------------------------------------------
function out = Vxx(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_atan(X,Y,Z,R);
    out = -G*rho * sum(S .* vals);
end

function out = Vyy(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_atan(Y,X,Z,R);
    out = -G*rho * sum(S .* vals);
end

function out = Vzz(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_atan(Z,X,Y,R);
    out = -G*rho * sum(S .* vals);
end


%% -----------------------------------------------
% Off-diagonal components: Vxy, Vxz, Vyz
%% -----------------------------------------------
function out = Vxy(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_ln(Z,R);
    out = -G*rho * sum(S .* (-vals));
end

function out = Vxz(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_ln(Y,R);
    out = -G*rho * sum(S .* (-vals));
end

function out = Vyz(x,y,z,L,B,D,rho,G)
    [X,Y,Z,S] = GravCuboidTensor.vertex_coordinates(x,y,z,L,B,D);
    R = sqrt(X.*X + Y.*Y + Z.*Z);
    vals = GravCuboidTensor.tensor_ln(X,R);
    out = -G*rho * sum(S .* (-vals));
end


%% -----------------------------------------------
% Full tensor (3×3)
%% -----------------------------------------------
function T = tensor_point(x,y,z,L,B,D,rho,G)
    T = [ ...
        GravCuboidTensor.Vxx(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vxy(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vxz(x,y,z,L,B,D,rho,G);
        GravCuboidTensor.Vxy(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vyy(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vyz(x,y,z,L,B,D,rho,G);
        GravCuboidTensor.Vxz(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vyz(x,y,z,L,B,D,rho,G), ...
        GravCuboidTensor.Vzz(x,y,z,L,B,D,rho,G) ];
end


%% -----------------------------------------------
% Batch tensor evaluation
%% -----------------------------------------------
function T = tensor_batch(points,L,B,D,rho,G,threshold)
    if nargin < 7
        threshold = 100;
    end

    pts = double(points);
    N = size(pts,1);

    T = zeros(N,3,3);

    if N > threshold
        parfor i = 1:N
            T(i,:,:) = GravCuboidTensor.tensor_point( ...
                    pts(i,1),pts(i,2),pts(i,3), L,B,D,rho,G);
        end
    else
        for i = 1:N
            T(i,:,:) = GravCuboidTensor.tensor_point( ...
                    pts(i,1),pts(i,2),pts(i,3), L,B,D,rho,G);
        end
    end
end

end % methods
end % classdef
