function write_frames_into_video(frames, vdofn)
%vdofn=[figsvdir_loc filesep VideoName '.mp4'];
Obj=VideoWriter(vdofn,'MPEG-4');
Obj.FrameRate=2;
Obj.Quality=100;
open(Obj)

for it = 1:length(frames)
    writeVideo(Obj, frames(it));
end
close(Obj)


return