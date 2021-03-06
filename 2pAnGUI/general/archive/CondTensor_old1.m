function [y1 y2] = CondTensor(Tlim,b,shiftflag)

%old1 does not have predelay/postdelay  or the image shift stuff

%Compute the tensor for each condition
%
%b is a 2D vector corresponding the the beginning and end of
%the baseline subtraction images, in milliseconds. e.g. varargin = {[0 500]} sums
%the images from 0 to .5 seconds for each repetition and then subtracts it
%from the mean response in the repeat.
%
%Tlim is like b, but for the range over which images are averaged

global ACQinfo bsflag bcond 

nc = pepgetnoconditions;
nr = pepgetnorepeats;

%Find blank condition:
bcond = NaN; 
for(i=0:nc-1)
    pepsetcondition(i)
    if(pepblank)       %Identify blank
       bcond = i;
       break
    end
end

imsize = ACQinfo.pixelsPerLine*ACQinfo.linesPerFrame;

pepsetcondition(0)

y1 = cell(1,nc);
y2 = cell(1,nc);

%Get sample period (ms/pixel)
sp = ACQinfo.msPerLine/ACQinfo.pixelsPerLine; %(msPerLine is actually sec/line)

for c = 0:nc-1
 
    if c == bcond

        y1{c+1} = [];
        y2{c+1} = [];

    else

        y1{c+1} = 0;
        y2{c+1} = 0;

        pepsetcondition(c);
        nr = pepgetnorepeats;

        for r = 0:nr-1
            pepsetcondition(c)
            pepsetrepeat(r);

            CHs = GetTrialData([1 1 1]);
            
            for z = 1:length(CHs{3}(1,1,:))
                synctcourse(:,:,z) = CHs{3}(:,:,z)';
            end
            synctimes = reconstructSync(synctcourse(:));

            idstart = round(synctimes(1)/sp) + 1;  %index of the first sync (should be ~-50)
            
            if bsflag == 1
                
                ttag = pepgettimetag;
                
                if ttag
                    pepsettimetag(ttag-1);
                    CHsblank = GetTrialData([1 1 0]);
                    
                    imstop = length(CHsblank{1}(1,1,:));
                    imstop = imstop-1;
                    imstart = imstop-6;

                    bimg1 = mean(CHsblank{1}(:,:,imstart:imstop),3);
                    bimg2 = mean(CHsblank{2}(:,:,imstart:imstop),3);
                else
                    pepsettimetag(ttag);
                
                    CHsblank = GetTrialData([1 1 0]);
                    imstop = length(CHsblank{1}(1,1,:))
                    imstart = 1;
                    imstop = 1;

                    bimg1 = mean(CHsblank{1}(:,:,imstart:imstop),3);
                    bimg2 = mean(CHsblank{2}(:,:,imstart:imstop),3);
                end
                
                for z = 1:length(CHs{1}(1,1,:))
                    CHs{1}(:,:,z) = CHs{1}(:,:,z) - bimg1;   %% baseline subtraction
                    CHs{2}(:,:,z) = CHs{2}(:,:,z) - bimg2;
                end
                
%                 for z = 1:length(CHs{1}(1,1,:))
%                     CHs{1}(:,:,z) = CHs{1}(:,:,z)./bimg1;   %% baseline subtraction
%                     CHs{2}(:,:,z) = CHs{2}(:,:,z)./bimg2;
%                 end

            end
            
            imstart = 1;
            imstop = 7;

            y1{c+1} = y1{c+1} + CHs{1}(:,:,imstart:imstop); %Add repeats
            y2{c+1} = y2{c+1} + CHs{2}(:,:,imstart:imstop); 
            clear CHs

        end
    end
end
