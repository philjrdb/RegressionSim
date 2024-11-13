function [c_idx] = consec_idx(indices,threshold)

% Derives logical index of consecutive indices >= threshold
% e.g. if indices = [1 3 5 6 7 8 10 11 12 15], threshold = 3
%           c_idx = [0 0 1 1 1 1  1  1  1  0]
%     i.e. indices(c_idx) = [5 6 7 8 10 11 12]

%  Copyright 2019 Philip Jean-Richard-dit-Bressel, UNSW Sydney

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

k = [true,diff(indices)~=1];
s = cumsum(k);
x =  histc(s,1:s(end));
idx = find(k);
consecutive = idx(x>=threshold);

c_idx = false(size(indices));

for c=1:length(consecutive)
   x_idx = consecutive(c) == idx;
   c_idx(consecutive(c):(consecutive(c)+x(x_idx)-1)) = true;
end  