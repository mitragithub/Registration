function frame2GifResample(frame,filename,time,desiredtime,flag,delaytime)
% if flag is 1, then go forwards and backwards
if nargin == 4
    flag = 0;
end
if flag
    frame = [frame(1:end),frame(end-1:-1:2)];
end
if nargin < 6
    delaytime = 0;
end


% keyboard
% make new frames, then call usual function
newframe = repmat(frame(1),[1,length(desiredtime)]);
for i = 1 : length(desiredtime)
    % find index
    ilow = find(desiredtime(i) >= time, 1, 'last');
    ihigh = ilow + 1;
    if ilow < 1
        ilow = 1;
    end
    if ilow > length(time)
        ilow = length(time);
    end
    if ihigh < 1
        ihigh = 1;
    end
    if ihigh > length(time)
        ihigh = length(time);
    end
    p = (desiredtime(i) - time(ilow))/(time(ihigh)-time(ilow));
    if isnan(p)
        p = 0;
    end
    if p == Inf
        p = 0;
    end
    newframe(i).cdata = uint8(double(frame(ilow).cdata)*(1-p) + double(frame(ihigh).cdata)*p);
end

 frame2Gif(newframe,filename,flag,delaytime);