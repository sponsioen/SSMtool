function C = sym2cell(S)
%SYM2CELL Convert sym array to cell array
%   C = sym2cell(S) converts a sym array S into cell array C by placing
%   each element of S into a separate cell in C. The output array has the
%   same size and dimensions as the input array.
%
%   See also CELL2SYM, CELL2MAT, MAT2CELL, NUM2CELL.

% Copyright 2015 The MathWorks, Inc.

C = num2cell(S);
