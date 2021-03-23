function []=plot_real_abs_fit(n_COM_alllaps,n_slope,pvalue,yrange,textpos)
sig_f_n_COM_alllaps = n_COM_alllaps(pvalue<0.05 & n_slope>0);
sig_f_n_slope = n_slope(pvalue<0.05 & n_slope>0);

sig_b_n_COM_alllaps = n_COM_alllaps(pvalue<0.05 & n_slope<0);
sig_b_n_slope = n_slope(pvalue<0.05 & n_slope<0);


figure;
%subplot(1,2,1)
hold on
plot(sig_f_n_COM_alllaps,sig_f_n_slope,'o')
plot(sig_b_n_COM_alllaps,sig_b_n_slope,'o')
title('novel')
xlabel('PF position')
ylabel('slope')

model_f_n=fitlm(sig_f_n_COM_alllaps,sig_f_n_slope)
plot(sig_f_n_COM_alllaps,model_f_n.('Coefficients').('Estimate')(2).*sig_f_n_COM_alllaps+ model_f_n.('Coefficients').('Estimate')(1));
%title(['slope = ' num2str(model_n.('Coefficients').('Estimate')(2)*6)])
f_summary=anova(model_f_n,'summary');
p_value=f_summary.('pValue')(2);
text(textpos,yrange(2)-1.8,['intercept ' num2str(model_f_n.('Coefficients').('Estimate')(1)) ' slope ' num2str(model_f_n.('Coefficients').('Estimate')(2))])
text(textpos,yrange(2)-1.2,['p-value' num2str(p_value)])
text(textpos,yrange(2)-0.6,['R-square ' num2str(model_f_n.Rsquared.Ordinary)]);
title ('real value slope vs position')

model_b_n=fitlm(sig_b_n_COM_alllaps,sig_b_n_slope)
plot(sig_b_n_COM_alllaps,model_b_n.('Coefficients').('Estimate')(2).*sig_b_n_COM_alllaps+ model_b_n.('Coefficients').('Estimate')(1));
%title(['slope = ' num2str(model_n.('Coefficients').('Estimate')(2)*6)])
b_summary=anova(model_b_n,'summary');
p_value=b_summary.('pValue')(2);
text(textpos,yrange(1)+1.8,['intercept ' num2str(model_b_n.('Coefficients').('Estimate')(1)) ' slope ' num2str(model_b_n.('Coefficients').('Estimate')(2))])
text(textpos,yrange(1)+1.2,['p-value' num2str(p_value)])
text(textpos,yrange(1)+0.6,['R-square ' num2str(model_b_n.Rsquared.Ordinary)]);
title ('real value slope vs position')

ylim(yrange)

% subplot(1,2,2)
% hold on;
% plot(n_COM_alllaps,abs(n_slope),'o')
% xlabel('PF position')
% ylabel('slope')
% model_n=fitlm(n_COM_alllaps,abs(n_slope));
% plot(n_COM_alllaps,model_n.('Coefficients').('Estimate')(2).*n_COM_alllaps+ model_n.('Coefficients').('Estimate')(1));
% %title(['slope = ' num2str(model_f.('Coefficients').('Estimate')(2)*6)])
% summary=anova(model_n,'summary');
% p_value=summary.('pValue')(2);
% text(textpos,yrange(1)+0.3,['intercept ' num2str(model_n.('Coefficients').('Estimate')(1)) ' slope ' num2str(model_n.('Coefficients').('Estimate')(2))])
% text(textpos,yrange(1)+0.2,['p-value' num2str(p_value)])
% text(textpos,yrange(1)+0.1,['R-square ' num2str(model_n.Rsquared.Ordinary)]);
% % text(textpos,yrange(1)+0.3,['intercept ' num2str(model_n.('Coefficients').('Estimate')(1)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(1)) ])
% % text(textpos,yrange(1)+0.2,['slope ' num2str(model_n.('Coefficients').('Estimate')(2)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(2)) ])
% % text(textpos,yrange(1)+0.1,['R-square ' num2str(model_n.Rsquared.Ordinary)]);title('F abs value slope vs posotopn')
% ylim(yrange)
% title ('abs value slope vs position')
% set(gcf, 'Position',  [100, 100, 1200, 500])





end
