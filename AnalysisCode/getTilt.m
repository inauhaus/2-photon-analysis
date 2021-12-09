function tilt = getTilt(anim,ExpLogArray)

dim = size(ExpLogArray);

for i = 2:dim(1)
   if strcmp(ExpLogArray{i,1},anim)
       tilt = str2num(char(ExpLogArray{i,34}));
       break
   end
    
end