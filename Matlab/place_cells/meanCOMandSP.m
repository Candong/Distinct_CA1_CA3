function [meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
binM=1:numbins;
COM=sum(binmean(:,start_b:end_b).*binM(start_b:end_b),2)./sum(binmean(:,start_b:end_b),2);
COM(isnan(COM))=0;
A=max(binmean(:,start_b:end_b),[],2);
meanCOM=sum(A.*COM)/sum(A);

SP=1/(sqrt(sum(A.*((COM-meanCOM).^2))/sum(A)));
 end 