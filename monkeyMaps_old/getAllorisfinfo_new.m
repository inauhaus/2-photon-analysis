global G_RChandles G_handles PW NB TC

%% ab2; u000_014

clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair

global maskS cellS

anim = 'ab2';
expt = 'u000_014';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

%The first experiment analyzed initializes these structures
%Get pairwise info
for i = 1:length(PW.dori)  %loop through each distance
    dori{1}{i} = PW.dori{i}(:);
    dsf{1}{i} = PW.dsf{i}(:);
    doriEuc{1}{i} = PW.doriEuc{i}(:);
    dsfEuc{1}{i} = PW.dsfEuc{i}(:);
    ax{1}{i} = PW.ax{i}(:);
    dist{1}{i} = PW.dist{i}(:);
    doripair{1}{i} = PW.doripair{i};
    dsfpair{1}{i} = PW.dsfpair{i};
    dtcoripair{1}{i} = PW.dtcoripair{i};
    dtcsfpair{1}{i} = PW.dtcsfpair{i};
end

%Initialize neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    end
end

%Initialize single-cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    else
        eval([fnames{i} '= dum;']);
    end
end

close all
%% ab2; u000_068

global maskS cellS

anim = 'ab2';
expt = 'u000_068';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{2}{i} = PW.dori{i}(:);
    dsf{2}{i} = PW.dsf{i}(:);
    doriEuc{2}{i} = PW.doriEuc{i}(:);
    dsfEuc{2}{i} = PW.dsfEuc{i}(:);
    ax{2}{i} = PW.ax{i}(:);
    dist{2}{i} = PW.dist{i}(:);
    doripair{2}{i} = PW.doripair{i};
    dsfpair{2}{i} = PW.dsfpair{i};
    dtcoripair{2}{i} = PW.dtcoripair{i};
    dtcsfpair{2}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Accumulate single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end

close all
%% ab2; u000_087

global maskS cellS

anim = 'ab2';
expt = 'u000_087';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[14:30]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{3}{i} = PW.dori{i}(:);
    dsf{3}{i} = PW.dsf{i}(:);
    doriEuc{3}{i} = PW.doriEuc{i}(:);
    dsfEuc{3}{i} = PW.dsfEuc{i}(:);
    ax{3}{i} = PW.ax{i}(:);
    dist{3}{i} = PW.dist{i}(:);
    doripair{3}{i} = PW.doripair{i};
    dsfpair{3}{i} = PW.dsfpair{i};
    dtcoripair{3}{i} = PW.dtcoripair{i};
    dtcsfpair{3}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:)]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end

end

close all
%% ab2; u001_004

global maskS cellS

anim = 'ab2';
expt = 'u001_004';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{4}{i} = PW.dori{i}(:);
    dsf{4}{i} = PW.dsf{i}(:);
    doriEuc{4}{i} = PW.doriEuc{i}(:);
    dsfEuc{4}{i} = PW.dsfEuc{i}(:);
    ax{4}{i} = PW.ax{i}(:);
    dist{4}{i} = PW.dist{i}(:);
    doripair{4}{i} = PW.doripair{i};
    dsfpair{4}{i} = PW.dsfpair{i};
    dtcoripair{4}{i} = PW.dtcoripair{i};
    dtcsfpair{4}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:)]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end

close all
%% ab3; u002_014

global maskS cellS

anim = 'ab3';
expt = 'u002_014';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{5}{i} = PW.dori{i}(:);
    dsf{5}{i} = PW.dsf{i}(:);
    doriEuc{5}{i} = PW.doriEuc{i}(:);
    dsfEuc{5}{i} = PW.dsfEuc{i}(:);
    ax{5}{i} = PW.ax{i}(:);
    dist{5}{i} = PW.dist{i}(:);
    doripair{5}{i} = PW.doripair{i};
    dsfpair{5}{i} = PW.dsfpair{i};
    dtcoripair{5}{i} = PW.dtcoripair{i};
    dtcsfpair{5}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:)]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end

close all
%% ab3; u002_041

global maskS cellS

anim = 'ab3';
expt = 'u002_041';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{6}{i} = PW.dori{i}(:);
    dsf{6}{i} = PW.dsf{i}(:);
    doriEuc{6}{i} = PW.doriEuc{i}(:);
    dsfEuc{6}{i} = PW.dsfEuc{i}(:);
    ax{6}{i} = PW.ax{i}(:);
    dist{6}{i} = PW.dist{i}(:);
    doripair{6}{i} = PW.doripair{i};
    dsfpair{6}{i} = PW.dsfpair{i};
    dtcoripair{7}{i} = PW.dtcoripair{i};
    dtcsfpair{7}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end
close all
%% aa9; u001_027

global maskS cellS

anim = 'aa9';
expt = 'u001_027';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{7}{i} = PW.dori{i}(:);
    dsf{7}{i} = PW.dsf{i}(:);
    doriEuc{7}{i} = PW.doriEuc{i}(:);
    dsfEuc{7}{i} = PW.dsfEuc{i}(:);
    ax{7}{i} = PW.ax{i}(:);
    dist{7}{i} = PW.dist{i}(:);
    doripair{7}{i} = PW.doripair{i};
    dsfpair{7}{i} = PW.dsfpair{i};
    dtcoripair{7}{i} = PW.dtcoripair{i};
    dtcsfpair{7}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end

close all

%% ab4; u001_003

global maskS cellS

anim = 'ab4';
expt = 'u001_003';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{8}{i} = PW.dori{i}(:);
    dsf{8}{i} = PW.dsf{i}(:);
    doriEuc{8}{i} = PW.doriEuc{i}(:);
    dsfEuc{8}{i} = PW.dsfEuc{i}(:);
    ax{8}{i} = PW.ax{i}(:);
    dist{8}{i} = PW.dist{i}(:);    
    doripair{8}{i} = PW.doripair{i};
    dsfpair{8}{i} = PW.dsfpair{i};
    dtcoripair{8}{i} = PW.dtcoripair{i};
    dtcsfpair{8}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end

close all
%% ab3; u002_027  (Why did I exclude him before?)
 
global maskS cellS

anim = 'ab3';
expt = 'u002_027';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Gkernelplots4

for i = 1:length(PW.dori)  %loop through each distance
    dori{9}{i} = PW.dori{i}(:);
    dsf{9}{i} = PW.dsf{i}(:);
    doriEuc{9}{i} = PW.doriEuc{i}(:);
    dsfEuc{9}{i} = PW.dsfEuc{i}(:);
    ax{9}{i} = PW.ax{i}(:);
    dist{9}{i} = PW.dist{i}(:);
    doripair{9}{i} = PW.doripair{i};
    dsfpair{9}{i} = PW.dsfpair{i};
    dtcoripair{9}{i} = PW.dtcoripair{i};
    dtcsfpair{9}{i} = PW.dtcsfpair{i};
end

%Accumulate neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    end
end

%Get single cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
    else
        eval([fnames{i} '= [' fnames{i} '; dum;]']);
    end
end
close all

%% ab3; u004_016  (Why did I exclude him before?)
 
% global maskS cellS
% 
% anim = 'ab3';
% expt = 'u004_016';
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% 
% maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
% maskpath = [maskroot anim '_' expt];
% traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
% tracepath = [traceroot anim '_' expt '_cellS'];
% 
% load(maskpath,'maskS') 
% load(tracepath,'cellS') 
% 
% set(G_RChandles.dropTrials,'string','[]')
% 
% Gkernelplots4
% 
% for i = 1:length(PW.dori)  %loop through each distance
%     dori{10}{i} = PW.dori{i}(:);
%     dsf{10}{i} = PW.dsf{i}(:);
%     doriEuc{10}{i} = PW.doriEuc{i}(:);
%     dsfEuc{10}{i} = PW.dsfEuc{i}(:);
%     ax{10}{i} = PW.ax{i}(:);
%     dist{10}{i} = PW.dist{i}(:);
%     doripair{10}{i} = PW.doripair{i};
%     dsfpair{10}{i} = PW.dsfpair{i};
% end
% 
% %Accumulate neighborhood info
% fnames = fieldnames(NB);
% for i = 1:length(fnames)
%     dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
%     if min(dim) == 1
%         eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
%     end
% end
% 
% %Get single cell info
% fnames = fieldnames(TC);
% for i = 1:length(fnames)
%     dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
%     if min(dim) == 1
%         eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
%     else
%         eval([fnames{i} '= [' fnames{i} '; dum;]']);
%     end
% end
% close all
% 
