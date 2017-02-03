function [new_points,N_matrix]=Normalize2Ddata(points)
% Isometric Normalization of 2D points
%
% Given a set of points in R^2, outputs a normalization matrix that, applied
% to the points (in homogeneous coordinates), transforms them into having 
% mean (0,0) and mean distance to the center equal to sqrt(2).
%
% Input arguments:
% points - 2xn - vector of n points of dimension 2
%
% Output arguments:
% N - 3x3 - isometric normalization matrix
% new_points - 3xn - vector of the n normalized points of dimension 2 in
%                    homogeneous coordinates

n=size(points,2);
points0=mean(points,2);
norm0=mean(sqrt(sum((points-repmat(points0,1,n)).^2,1)));
N_matrix=diag([sqrt(2)/norm0;sqrt(2)/norm0;1]);
N_matrix(1:2,3)=-sqrt(2)*points0/norm0;

new_points=N_matrix(1:2,:)*[points;ones(1,n)];