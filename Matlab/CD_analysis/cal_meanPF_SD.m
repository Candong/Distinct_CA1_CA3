function SD=cal_meanPF_SD(mean_binmean)
% in cm
xi=1:length(mean_binmean);
xi=xi'*6;
COM=sum(xi.*(mean_binmean/sum(mean_binmean)));
SD=sqrt(sum((xi-COM).^2.*(mean_binmean/sum(mean_binmean))));
%skew = sum(((xi-COM)/sig).^3.*(mean_binmean/sum(mean_binmean)));

end