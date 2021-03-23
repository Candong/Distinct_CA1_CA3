function result=isclipped(binmean)
% find out if the PF is clipped.
activebin=find(mean(binmean,2)>0);
startbin=activebin(1);
endbin=activebin(end);
[meanCOM SP]=meanCOMandSP(binmean',1,size(binmean,1),size(binmean,1));

if startbin==1 & (endbin-meanCOM)>(meanCOM-startbin)
    result=1;
elseif endbin==size(binmean,1) & (endbin-meanCOM)<(meanCOM-startbin)
    result = 1;
else
    result=0;
end


end