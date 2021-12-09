function CHs = GetTrialData(chvec,tcr)

%example: GetTrialData([1 1 0],[cond rep]); or GetTrialData([1 1 0],trial)

%if tcr has one element, then it is the trial no.  if it has two
%elements then it is the cond/repeat

global twophDATADIR AUE ACQinfo

if length(tcr) == 1
    trial = tcr;
elseif length(tcr) == 2
    cond = tcr(1);
    rep = tcr(2);
    trial = gettrial(cond,rep);
end

filepath = [twophDATADIR AUE ' ' sprintf('%03d',trial) '.tif'];

tf = imformats('tif');
info = feval(tf.info, filepath);

Nimages = length(info);

framestart = 1;
framestop = Nimages;

if ACQinfo.acquiringChannel3
    dFrame = 3;
else
    dFrame = 2;
end

dim = [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine];

i = 1;
if chvec(1) == 1
    k = 1;
    dim2 = [dim length(framestart:dFrame:framestop)];
    CHs{i} = zeros(dim2,'single');
    for frame=framestart:dFrame:framestop
        A = feval(tf.read, filepath, frame);
        CHs{i}(:,:,k) = single(A);
        k = k+1;
    end
    i = 2;
end

if chvec(2) == 1
    k = 1;
    dim2 = [dim length((framestart+1):dFrame:framestop)];
    CHs{i} = zeros(dim2,'single');
    for frame=(framestart+1):dFrame:framestop
        A = feval(tf.read, filepath, frame);
        CHs{i}(:,:,k) = single(A);
        k = k+1;
    end
    i = i+1;
end

if chvec(3) == 1
    k = 1;    
    dim2 = [dim length((framestart+2):dFrame:framestop)];
    CHs{i} = zeros(dim2,'single');    
    for frame=(framestart+2):dFrame:framestop
        A = feval(tf.read, filepath, frame);
        CHs{i}(:,:,k) = single(A);
        k = k+1;
    end
end

clear A

