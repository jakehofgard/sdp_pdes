function [imon] = imom(ocp, mon)
%@POCP/PRIVATE/IMON - Internal use only.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 by Carlo Savorgnan 
% 
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License version 2 as 
% published by the Free Software Foundation. 
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License 
% along with this program; if not, write to the 
% Free Software Foundation, Inc., 
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Substitutes in the vector of monomials mon the variables.
%   E.g. suppose mom is the monomial x(1)^2*x(2) where x is the state
%   vector. Furthermore, assume that the value of x(1) at time zero has
%   been assigned, while x(2) has to be determined. In the returned fmon
%   x(1)^2 will be substituted with the numeric value of the second order
%   moment of the variable. x(2) will be substituted with the corresponding
%   variable in the occupation measure at time 0.

imon = mon;

for mind = 1:length(imon)
    list = listvar(imon(mind));
    for vind = 1:length(list)
        ex = deg(imon(mind), list(vind));
        if isequal(list(vind), ocp.time)
            %time variable
            imon(mind) = subs(imon(mind), list(vind), 0);
        elseif ~isempty(locate(ocp.icon.dirac.var, list(vind)))
            %Dirac's delta
            cind = locate(ocp.icon.dirac.var, list(vind));
            tmp = 0;
            for pind = 1:length(ocp.icon.dirac.weight)
                tmp = tmp + ocp.icon.dirac.weight(pind)*ocp.icon.dirac.value(cind, pind)^ex;
            end
            imon(mind) = subs(imon(mind), list(vind), 1)*tmp;
        elseif ~isempty(locate(ocp.icon.unif.var, list(vind)))
            %uniform distribution
            cind = locate(ocp.icon.unif.var, list(vind));
            tmp = 1/(ex+1)*(ocp.icon.unif.interval(cind, 2)^(ex+1) - ocp.icon.unif.interval(cind, 1)^(ex+1))/(ocp.icon.unif.interval(cind, 2)-ocp.icon.unif.interval(cind, 1));
            imon(mind) = subs(imon(mind), list(vind), 1)*tmp;
        else
            %the statistics of the variable hasn't been assigned
            imon(mind) = subs(imon(mind), list(vind), ocp.istate(locate(ocp.state, list(vind))));
        end
    end
end
            