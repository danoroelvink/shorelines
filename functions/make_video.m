function make_video(S,vi)
% function make_video(S,vi)
%
% This function writes video frames
%
% INPUT:
%    vi   input fgure frame used for the video
%    S
%         .video        : switch for making videos (0/1)
%         .outputdir    : name of the output directory (as string)
%         .diffraction  : switch for using diffraction (alternatively used for naming the output directory)
%         .xhard        : x-coordinates of the coastal structures (alternatively used for naming the output directory)
%         .yhard        : y-coordinates of the coastal structures (alternatively used for naming the output directory)
%         .phiw0        : offshore wave direction, static parameter (alternatively used for naming the output directory)
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    if S.video
        % label=strcat(S.trform,'_',S.boundaryconditionstart,'_ds0_',num2str(S.ds0),'_tc_',num2str(S.spread),'_phiw_',num2str(S.phiw0));
        if  S.diffraction==1
            label=strcat('xhrd_',num2str(S.xhard),'_yhrd_',num2str(S.yhard),'_phi0_',num2str(S.phiw0));
            % label=strcat('spit_',num2str(S.spitwidth),'ds0_',num2str(S.ds0),'_phi0_',num2str(S.phiw0));
        else
            label=[S.outputdir,'\animation'];
        end
        FR=10;
        % video=VideoWriter(label,'MPEG-4');
        if ~exist(S.outputdir,'dir')
            mkdir(S.outputdir);
        end
        
        video=VideoWriter(label,'Motion JPEG AVI'); % for surfsara runs
        video.FrameRate=FR;
        open(video)
        try
            writeVideo(video,vi)
        end
        close(video)
    end
end    
