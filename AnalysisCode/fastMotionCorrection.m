function CH = fastMotionCorrection(CH)

global cellS

%% Coarse x/y aligmnment using the motion info in cellS. Uses template match.

mbest = cellS.motionInfo.mbest{t}(1:CHdim(3));
nbest = cellS.motionInfo.nbest{t}(1:CHdim(3));
CH = align2Template(CH,mbest,nbest);

%% Shift correction using optic flow. Refinement.

%sigx = [3 1];
sigx = 2;
[CH dx dy] = TensorOpticFlowRigidCorrection(CH, sigx);


%% Get motion information and subtract it [using the unsmoothed version]

[CH dx dy] = TensorSubtractFlow(CH);
