function M=crossM(v)
% CROSSM(v) produces the 3x3 matrix corresponding to the cross product of
% vector v.
%
% Copyright (c) 2017 Laura F. Julia
M=[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];