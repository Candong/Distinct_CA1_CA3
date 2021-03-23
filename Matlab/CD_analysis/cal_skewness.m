function skew=cal_skewness(mean_binmean)

xi=1:length(mean_binmean);
xi=xi';
COM=sum(xi.*(mean_binmean/sum(mean_binmean)));
sig=sqrt(sum((xi-COM).^2.*(mean_binmean/sum(mean_binmean))));
skew = sum(((xi-COM)/sig).^3.*(mean_binmean/sum(mean_binmean)));

end