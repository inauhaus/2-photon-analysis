pF0

global G_handles Analyzer cellS maskS idExamp alignTemp ACQinfo

set(G_handles.slowMotionFlag,'Value',1); 
set(G_handles.searchRange,'String','9'); %movement correction search range 

%I have data spread across multiple drives. It will look in each one...
dataRootA = 'g:\2p_data\';
dataRootB = 'h:\2p_data\';
dataRootC = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%% 
clear animAll expt_mask exvec expts

idx = 1;

animAll{idx} = 'rl8'; %1/22/18
expt_mask{idx} = 'u001_002'; %ori experiment
exvec{idx} = [11 9   8   7   6   5   [4 10]];

idx = idx+1;
animAll{idx} = 'rl9'; %1/22/18
expt_mask{idx} = 'u002_003'; %ori experiment
exvec{idx} = [2 10 9 7 6 5 4 11];

idx = idx+1;
animAll{idx} = 'rl7'; %1/22/18
expt_mask{idx} = 'u001_002'; %ori experiment
exvec{idx} = [1 8 7 6 5 4 3];

idx = idx+1;
animAll{idx} = 'rl7'; %1/22/18
expt_mask{idx} = 'u001_013'; %ori experiment
exvec{idx} = [12 19 18 17 16 15 14];
 
idx = idx+1;
animAll{idx} = 'rl5'; %1/22/18
expt_mask{idx} = 'u001_003'; %ori experiment
exvec{idx} = [2 11 10 9 8 7 4 5 12];
 
idx = idx+1;
animAll{idx} = 'rl7'; %1/22/18
expt_mask{idx} = 'u002_003'; %ori experiment
exvec{idx} = [21 16 17 14 15 12 13 10 11 9 18 5 8 20 6 7 22] ;

idx = idx+1;
animAll{idx} = 'rm2'; %1/22/18
expt_mask{idx} = 'u000_002'; %ori experiment
exvec{idx} = [1 13  11  10  9  8  7   6   [4 5]];



idx = idx+1;
animAll{idx} = 'rm3'; %1/22/18
expt_mask{idx} = 'u000_002'; %ori experiment
exvec{idx} = [1 [17 18] [15 16] 14  12  [10 11] [8 9]   [6 7]   [3 4]];

idx = idx+1;
animAll{idx} = 'rm4'; %1/22/18
expt_mask{idx} = 'u001_006'; %ori experiment
exvec{idx} = [5     15  14  13  12  11  [9 10]  [7 8]];

idx = idx+1;
animAll{idx} = 'rm6'; %1/22/18
expt_mask{idx} = 'u001_003'; %ori experiment
exvec{idx} = [1 12  11  10  9   8   7   6   5];

% idx = idx+1;
% animAll{idx} = 'rm3'; %1/22/18
% expt_mask{idx} = 'u001_004'; %Mine
%%expt_mask{idx} = 'u001_003'; %Issac's 
% exvec{idx} = [2   12  11  10  9   8   [6 7]   5   4];

% idx = idx+1;
% animAll{idx} = 'rm4'; %1/22/18
% expt_mask{idx} = 'u002_005';  %Mine 
%%expt_mask{idx} = 'u002_003'; %Issac's 
% exvec{idx} = [2   12  11  10  9   8   7   6   5];

idx = idx+1;
animAll{idx} = 'rm5'; %1/22/18
expt_mask{idx} = 'u002_002'; %ori experiment
exvec{idx} = [3 16  15  11  10  9   8   5   4];

idx = idx+1;
animAll{idx} = 'rn1'; %1/22/18
expt_mask{idx} = 'u001_003'; %ori experiment
exvec{idx} = [20    14  13  12  11  10  9   [7 8]   [5 19]];

idx = idx+1;
animAll{idx} = 'rn2'; %1/22/18
expt_mask{idx} = 'u000_002'; %ori experiment
exvec{idx} = [1 [17 18] [14 15] 13  12  [10 11] [8 9]   [6 7]   [3 4]];

idx = idx+1;
animAll{idx} = 'rn4'; %1/22/18
expt_mask{idx} = 'u001_007'; %ori experiment
exvec{idx} = [5 7 8 9 10 11 12 17];

%%
clear expts
for i = 1:length(animAll)
    i
   for j = 1:length(exvec{i})
      expts{i}{j} = ['u00' num2str(expt_mask{i}(4)) '_' sprintf('%03d',exvec{i}(j))]; 
   end
end


%%
exdom = 1:length(animAll);

%for e = 1:length(exdom)  
for e = 1:1  %Still need 10 12 13
    ex = exdom(e);
    
    anim = animAll{ex};
    
    %Load mask
    expt = expt_mask{ex};
    maskroot = 'C:\CellMasks_Issac\';
    %maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    load(maskpath,'maskS')
    
    %Make alignment template using mask experiment
    try
        dum = single(sbxread([dataRootA anim '\' anim '_' expt_mask{ex}(2:end)],1,100));
        dataRoot = dataRootA;
    catch
        try
            dum = single(sbxread([dataRootB anim '\' anim '_' expt_mask{ex}(2:end)],1,100));
            dataRoot = dataRootB;
        catch
            dum = single(sbxread([dataRootC anim '\' anim '_' expt_mask{ex}(2:end)],1,100));
            dataRoot = dataRootC;
        end
    end        
        
    'using first 100 frames'
    alignTemp = median(squeeze(dum),3);
    alignTemp = alignTemp(:,54:729);
    
    %% Loop through experiments: make traces and save
    for edum = 1:length(expts{ex})
        expt = expts{ex}{edum};
        set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])

        set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
        Gsetdirectories %Load raw experiment and analyzer file
        
        fixSBsyncs
        %%%%%%%%%%%
        
        %Make traces for exp x
        processButton(0)
        
        %Save traces for exp x
        traceroot = 'C:\CellTraces\';
        tracepath = [traceroot anim '_' expt '_cellS'];
        save([tracepath '_aligned to ' expt_mask{ex}(2:end)],'cellS')
        
    end
    
    
end







%% Old... 
% 
% 
% 
% 
% 
% pF0
% 
% global G_handles Analyzer cellS maskS idExamp alignTemp ACQinfo
% 
% set(G_handles.slowMotionFlag,'Value',1); %baseline subtraction
% set(G_handles.searchRange,'String','9'); %baseline subtraction
% 
% dataRoot = 'g:\2p_data\';
% anaRoot = 'c:\AnalyzerFiles\';
% 
% %%
% 
% %this experiment sucks
% animAll{1} = 'rk9';
% expt_mask{1} = 'u001_002';
% expt_0v90{1} = 'u001_007';  %could try 4 and 5
% expt_45v135{1} = 'u001_008'; %these suck
% expt_Kal{1} = 'u001_001';
% idExampAll{1} = [];
% 
% animAll{2} = 'rm0';
% expt_mask{2} = 'u000_007';
% expt_0v90{2} = 'u000_008';
% expt_45v135{2} = 'u000_007';
% expt_Kal{2} = 'u000_005';
% 
% %
% % animAll{1} = 'rl9';
% % expt_mask{1} = 'u001_004';
% % expt_0v90{1} = 'u001_008';
% % expt_45v135{1} = 'u001_007';
% % expt_Kal{1} = 'u001_005';
% 
% animAll{1} = 'rm2';
% expt_mask{1} = 'u001_008';
% expt_0v90{1} = 'u001_007';
% expt_45v135{1} = 'u001_008';
% expt_Kal{1} = 'u001_002';
% 
% animAll{1} = 'rm5';    %Where is he?
% expt_mask{1} = 'u001_007';
% expt_0v90{1} = 'u001_006';
% expt_45v135{1} = 'u001_007';
% expt_Kal{1} = 'u001_002';
% 
% %new new
% 
% animAll{1} = 'rm5';
% expt_mask{1} = 'u002_007';
% expt_0v90{1} = 'u002_006';
% expt_45v135{1} = 'u002_007';
% expt_Kal{1} = 'u002_003';
% 
% animAll{2} = 'rm2';
% expt_mask{2} = 'u002_005';
% expt_0v90{2} = 'u002_004';
% expt_45v135{2} = 'u002_005';
% expt_Kal{2} = 'u002_001';
% 
% animAll{3} = 'rr1';
% expt_mask{3} = 'u001_008';
% expt_0v90{3} = 'u001_007';
% expt_45v135{3} = 'u001_008';
% expt_Kal{3} = 'u001_001';
% 
% animAll{4} = 'rr0';
% expt_mask{4} = 'u001_008';
% expt_0v90{4} = 'u001_007';
% expt_45v135{4} = 'u001_008';
% expt_Kal{4} = 'u001_003';
% 
% animAll{5} = 'rr2';
% expt_mask{5} = 'u001_007';
% expt_0v90{5} = 'u001_006';
% expt_45v135{5} = 'u001_007';
% expt_Kal{5} = 'u001_002';
% 
% %%%new new
% 
% animAll{1} = 'rr9';
% expt_mask{1} = 'u001_005';
% expt_0v90{1} = 'u001_004';
% expt_45v135{1} = 'u001_005';
% expt_Kal{1} = 'u001_002';
% 
% animAll{2} = 'rr8';
% expt_mask{2} = 'u001_004';
% expt_0v90{2} = 'u001_003';
% expt_45v135{2} = 'u001_004';
% expt_Kal{2} = 'u001_001';
% 
% animAll{3} = 'rw5';
% expt_mask{3} = 'u001_009';
% expt_0v90{3} = 'u001_007';
% expt_45v135{3} = 'u001_009';
% expt_Kal{3} = 'u001_003';
% 
% animAll{1} = 'rw9';  
% expt_mask{1} = 'u001_006';  %looks like a lot of brain movement
% expt_0v90{1} = 'u001_005';
% expt_45v135{1} = 'u001_006';
% expt_Kal{1} = 'u001_002';
% 
% animAll{2} = 'rx0';
% expt_mask{2} = 'u001_005';
% expt_0v90{2} = 'u001_004';
% expt_45v135{2} = 'u001_005';
% expt_Kal{2} = 'u001_002';
% 
% animAll{3} = 'rw8';
% expt_mask{3} = 'u001_010';
% expt_0v90{3} = 'u001_009';
% expt_45v135{3} = 'u001_010';
% expt_Kal{3} = 'u001_007';
% 
% animAll{4} = 'rw5';
% expt_mask{4} = 'u002_005';
% expt_0v90{4} = 'u002_004';
% expt_45v135{4} = 'u002_005';
% expt_Kal{4} = 'u002_002';
% 
% animAll{5} = 'rx1';
% expt_mask{5} = 'u002_005';
% expt_0v90{5} = 'u002_004';
% expt_45v135{5} = 'u002_005';
% expt_Kal{5} = 'u002_002';
% 
% 'rm5_u002'
% 'rm2_u002'
% 'rr1_u001'
% 'rr0_u001'
% 'rr2_u001'
% 'rr9_u001'
% 'rr8_u001'
% 'rw5_u001'
% 'rw9_u001'
% 'rx0_u001'
% 'rw8_u001'
% 'rw5_u002'
% 'rx1_u002'
% 
% %%
% exdom = 1:1;
% 
% for e = 1:length(exdom)
%     
%     ex = exdom(e);
%     
%     anim = animAll{ex};
%     
%     %Load mask
%     expt = expt_mask{ex};
%     maskroot = 'C:\CellMasks\';
%     maskpath = [maskroot anim '_' expt(1:8)];
%     load(maskpath,'maskS')
%     
%     %Make alignment template using mask experiment
%     dum = single(sbxread([dataRoot anim '\' anim '_' expt_mask{ex}(2:end)],1,100));
%     'using first 100 frames'
%     alignTemp = median(squeeze(dum),3);
%     alignTemp = alignTemp(:,ACQinfo.unblanked);
%     
%     %% Load exp1, make traces, save
%     expt = expt_0v90{ex}; %S vs M
%     set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
%     set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
%     try
%         Gsetdirectories %Load raw experiment and analyzer file
%     catch
%         loadAnalyzer %Just load analyzer file
%     end
%     fixSBsyncs
%     %%%%%%%%%%%
%     
%     %Make traces for exp1
%     processButton(0)
%     
%     %Save traces for exp1
%     traceroot = 'C:\CellTraces\';
%     tracepath = [traceroot anim '_' expt '_cellS'];
%     save([tracepath '_aligned to ' expt_mask{ex}(2:end)],'cellS')
%     
%     %% Load exp2, make traces, save
%     expt = expt_45v135{ex}; %S vs M
%     set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
%     set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
%     try
%         Gsetdirectories %Load raw experiment and analyzer file
%     catch
%         loadAnalyzer %Just load analyzer file
%     end
%     fixSBsyncs
%     %%%%%%%%%%%
%     
%     %Make traces for exp1
%     processButton(0)
%     
%     %Save traces for exp2
%     traceroot = 'C:\CellTraces\';
%     tracepath = [traceroot anim '_' expt '_cellS'];
%     save([tracepath '_aligned to ' expt_mask{ex}(2:end)],'cellS')
%     
%     %% Load Kalatsky, make traces, save
%     expt = expt_Kal{ex}; %Kalatsky
%     set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
%     set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
%     try
%         Gsetdirectories %Load raw experiment and analyzer file
%     catch
%         loadAnalyzer %Just load analyzer file
%     end
%     fixSBsyncs
%     %%%%%%%%%%%
%     
%     %Make traces for exp1
%     processButton(0)
%     
%     %Save traces for exp2
%     traceroot = 'C:\CellTraces\';
%     tracepath = [traceroot anim '_' expt '_cellS'];
%     save([tracepath '_aligned to ' expt_mask{ex}(2:end)],'cellS')
%     
%     
% end
% 
% %%
% 'rk5_u001'*
% 'rk7_u001'*
% 'rk3_u002'*
% 'rl0_u001'*
% 'rl2_u000'*
% 'rl1_u002'*
% 'rk9_u001'  (Turrell)
% 'rl0_u002'*
% 'rk5_u002'*
% 'rm0_u000'   (Turrell)
% 'rl9_u001'  (Turrell)
% 'rm2_u001'  (Turrell)
% 'rm5_u001' (Turrell)
% 'rm5_u002'
% 'rm2_u002'
% 'rr1_u001'
% 'rr0_u001'
% 'rr2_u001'
% 'rr9_u001'
% 'rr8_u001'
% 'rw5_u001'
% 'rw9_u001'
% 'rx0_u001'
% 'rw8_u001'
% 'rw5_u002'
% 'rx1_u002'
% 
% %rk4 unit 1?
% %%
% animAll{1} = 'rk5';
% expt_mask{1} = 'u001_010';
% expt_0v90{1} = 'u001_010';
% expt_45v135{1} = 'u001_011';
% expt_Kal{1} = 'u001_002';
% idExampAll{1} = [3 21 75 81 92];
% 
% animAll{2} = 'rk7';
% expt_mask{2} = 'u001_006';
% expt_0v90{2} = 'u001_006';
% expt_45v135{2} = 'u001_007';
% expt_Kal{2} = 'u001_000';
% idExampAll{2} = [3 21];
% 
% animAll{3} = 'rl0';
% expt_mask{3} = 'u001_005';
% expt_0v90{3} = 'u001_005';
% expt_45v135{3} = 'u001_006';
% expt_Kal{3} = 'u001_001';
% idExampAll{3} = [19 21 25 49 58 78 98];
% 
% animAll{4} = 'rl1';
% expt_mask{4} = 'u002_006';
% expt_0v90{4} = 'u002_006';
% expt_45v135{4} = 'u002_007';
% expt_Kal{4} = 'u002_003';
% idExampAll{4} = [];
% 
% animAll{5} = 'rk3';
% expt_mask{5} = 'u002_009';
% expt_0v90{5} = 'u002_010';
% expt_45v135{5} = 'u002_011';
% expt_Kal{5} = 'u002_004';
% idExampAll{5} = [];
% 
% animAll{6} = 'rl2';
% expt_mask{6} = 'u000_005';
% expt_0v90{6} = 'u000_006';
% expt_45v135{6} = 'u000_007';
% expt_Kal{6} = 'u000_001';
% idExampAll{6} = [];
% 
% 
% animAll{7} = 'rk5';
% expt_mask{7} = 'u002_002';
% expt_0v90{7} = 'u002_004';
% expt_45v135{7} = 'u002_005';
% expt_Kal{7} = 'u002_001';
% idExampAll{7} = [];
% 
% 
% animAll{8} = 'rl0';
% expt_mask{8} = 'u002_006';
% expt_0v90{8} = 'u002_004';
% expt_45v135{8} = 'u002_006';
% expt_Kal{8} = 'u002_001';
% idExampAll{8} = [];
% 
% 
% animAll{9} = 'rk4';
% expt_mask{9} = 'u001_003'; 
% expt_0v90{9} = 'u001_005';
% expt_45v135{9} = 'u001_002';
% expt_Kal{9} = 'u001_002';
% idExampAll{9} = [];