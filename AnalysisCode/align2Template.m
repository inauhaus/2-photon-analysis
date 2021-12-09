function T = align2Template(T,mbest,nbest)

global G_handles cellS

fastMo = get(G_handles.fastXcorr,'Value');  %lateral movement correction

%frames = Fwin(1):Fwin(end);

if fastMo
    
    if isfield(cellS.motionInfo,'mbestFrame')
        
        for q = 1:size(T,3) %loop each frame
           % mbest = cellS.motionInfo.mbestFrame(frames(q));
            %nbest = cellS.motionInfo.nbestFrame(frames(q));
            
            %T(:,:,q) = circshift(T(:,:,q),[round(-mbest) round(-nbest) 0]);            
            T(:,:,q) = circshift_continous2(T(:,:,q),-nbest(q),-mbest(q));
            
        end
    else
        warning('Not applying motion correction.  Need to first "Get Motion Information".')
    end
end
