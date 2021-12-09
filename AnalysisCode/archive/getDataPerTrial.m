function CHs = getDataPerTrial(chvec,trial)
 
%Loads the data from a single trial (corresponds to 1 .tiff file) 

%This is different from getTrialData, in that you specify the trial number,
%instead of cond/repeat... you don't set condition and repeat before
%calling this function, but the experiment must be loaded.

global twophDATADIR pepLOADED pepANA

anim = pepANA.config.animal;

filepath = [twophDATADIR '\' anim '_' pepLOADED ' ' sprintf('%03d',trial) '.tif'];
%filepath = [twophDATADIR '\' pepLOADED ts '.tif'];

tf = imformats('tif');
info = feval(tf.info, filepath);

Nimages = length(info);

tf = imformats('tif');

framestart = 1;
framestop = Nimages;

i = 1;
if chvec(1) == 1
    k = 1;
    for frame=framestart:3:framestop
        A = feval(tf.read, filepath, frame);
        %A = A';
        %A = A(:);
        %CHs{i}(:,k) = double(A);
        CHs{i}(:,:,k) = double(A);
        k = k+1;
    end
    i = 2;
end

if chvec(2) == 1
    k = 1;
    for frame=(framestart+1):3:framestop
        A = feval(tf.read, filepath, frame);
        %A = A';
        %A = A(:);
        %CHs{i}(:,k) = double(A);
        CHs{i}(:,:,k) = double(A);
        k = k+1;
    end
    i = i+1;
end
if chvec(3) == 1
    k = 1;
    for frame=(framestart+2):3:framestop
        A = feval(tf.read, filepath, frame);
        %A = A';
        %A = A(:);
        %CHs{i}(:,k) = double(A);
        CHs{i}(:,:,k) = double(A);
        k = k+1;
    end
end

clear A

