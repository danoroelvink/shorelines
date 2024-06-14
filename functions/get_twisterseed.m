function [ret]=get_twisterseed(SEED)
% function [ret]=get_twisterseed(SEED)
% 
% INPUT:
%   SEED           seeding number used as basis for randomly generating the 'ret' 
%                  by using seem 'SEED' you nay reproduce the same sequence of random numbers
%
% OUTPUT:
%   ret            randomly generated number
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
%
%       Johan Reyns
%       j.reyns@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

   ret = uint32(zeros(625,1));
   ret(1) = SEED;
   for N = 1:623
       % initialize_generator
       % bit-xor (right shift by 30 bits)
       a=uint64(1812433253)*uint64(bitxor(ret(N),bitshift(ret(N),-30)))+N; % has to be uint64, otherwise in 4th iteration hit maximum of uint32!
       ret(N+1) = uint32(bitand(a,uint64(intmax('uint32')))); % untempered numbers
   end
   ret(end) = 1;   
end