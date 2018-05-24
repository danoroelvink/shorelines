function []=make_video(S,vi)
if S.video
% label=strcat(S.trform,'_',S.boundary_condition_start,'_ds0_',num2str(S.ds0),'_tc_',num2str(S.spread),'_phiw_',num2str(S.phiw0));
% label=strcat('xhrd_',num2str(S.x_hard),'_yhrd_',num2str(S.y_hard),'_phi0_',num2str(S.phiw0));
%label=strcat('spit_',num2str(S.spit_width),'ds0_',num2str(S.ds0),'_phi0_',num2str(S.phiw0));
label=[S.outputdir,'\animation'];
FR=10;
% video=VideoWriter(label,'MPEG-4');
video=VideoWriter(label,'Motion JPEG AVI'); % for surfsara runs
video.FrameRate=FR;
open(video)
writeVideo(video,vi)
close (video)
end