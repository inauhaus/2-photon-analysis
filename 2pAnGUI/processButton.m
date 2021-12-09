function processButton(F0im_flag)

global f0m f0m_var Tens Tens_var bsflag Flim repDom G_handles Analyzer

bsflag = 0;

t0 = cputime;

varflag = get(G_handles.EbarFlag,'Value');
bsflag = get(G_handles.basesub,'Value');
Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); %Frame stop in ms (to average)
b = str2double(get(G_handles.bstart,'String')); %in msec as well
b(2) = str2double(get(G_handles.bstop,'String'));

repdum = get(G_handles.repDom,'string');
if strcmp(repdum,'All')
    repDom = 1:getnorepeats(1);
else
    eval(['repDom = [' repdum '];'])
end

set(G_handles.status,'string','Processing...'), drawnow

if F0im_flag %Get F0 images (and possibly cell mask traces)
    [Tens Tens_var] = CondTensor5(b);  %%Compute entire space time block for each condition
    f0m = CondF0(Tens,Flim);   %%%Compute mean over time interval [Flim(1) Flim(2)]%%%
    %f0m_var = CondF0(Tens_var,Flim);
else %only get cell mask traces (for all trials)
      
    %CondMaskdata(slowMo,fastMo)
    CondMaskdata2

    getCellStats_NaN  %make 'cellS'
    
end


set(G_handles.status,'string','Done'), drawnow
%sound(.6*sin(2*pi*400/(1200)*(0:400)),1200)  %Signal done

t1 = cputime-t0;

set(G_handles.time,'string',num2str(t1))

set(G_handles.loaded,'string',[Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt])

set(G_handles.plot,'enable','on')

set(G_handles.save,'enable','on')


