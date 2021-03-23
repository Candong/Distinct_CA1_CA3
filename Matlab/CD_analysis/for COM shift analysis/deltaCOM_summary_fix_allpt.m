clear all; close all;

usecolor1=[0/255.0,113/255.0,188/255.0]; %blue
usecolor2='red'
usecolor3=[102/255.0, 45/255.0, 145/255.0]; %purple
usecolor4=[34/255.0,181/255.0,115/255.0]; % green
usecolor5=[241/255.0,90/255.0,36/255.0]; %orange

[filepaths, temp]=sort(uigetfile('*.mat', 'Chose deltaCOM files to load:','MultiSelect','on'));

for f=1:size(filepaths,2)
 if contains(filepaths{f},'_fnovel')
      if contains(filepaths{f},'CA1_')
        load(filepaths{f}); 
        CA1_f_x=f_x;
        CA1_nday1_x=n_x;
        CA1_f_y=f_y;
        CA1_nday1_y=n_y;
%         CA1_f_x_M=f_x_M;
%         CA1_f_y_M=f_y_M;
%         CA1_nday1_x_M=n_x_M;
%         CA1_nday1_x_M=n_y_M;
        
      elseif contains(filepaths{f},'CA3_')
        load(filepaths{f}); 
        CA3_f_x=f_x;
        CA3_nday1_x=n_x;
        CA3_f_y=f_y;
        CA3_nday1_y=n_y;
%         CA3_f_x_M=f_x_M;
%         CA3_f_y_M=f_y_M;
%         CA3_nday1_x_M=n_x_M;
%         CA3_nday1_x_M=n_y_M;        
%      
      end
      
 elseif contains(filepaths{f},'_novel2day')
     
      if contains(filepaths{f},'CA1_')
        load(filepaths{f}); 
        CA1_nday2_x=n_x;
        CA1_nday2_y=n_y;
%         CA1_nday2_x_M=n_x_M;
%         CA1_nday2_y_M=n_y_M;
      elseif contains(filepaths{f},'CA3_')
        load(filepaths{f}); 
        CA3_nday2_x=n_x;
        CA3_nday2_y=n_y;
%         CA3_nday2_x_M=n_x_M;
%         CA3_nday2_y_M=n_y_M;
        
     
      end
     
  
 end
end
%%
figure;
subplot(1,2,1)
[lap_range CA1_mean_fy CA1_f_sem]=cal_mean_sem(CA1_f_x, CA1_f_y);
hold on;
%plot(range,mean_y,'o')
errorbar(lap_range(1:end-1),CA1_mean_fy(1:end-1),CA1_f_sem(1:end-1),'o','linewidth',0.5,'color', usecolor1,'MarkerFaceColor',usecolor1)
%scatter(CA1_f_x, CA1_f_y,'filled','MarkerFaceAlpha',0.01)


[lap_range CA1_mean_nday1y CA1_nday1_sem]=cal_mean_sem(CA1_nday1_x, CA1_nday1_y);
errorbar(lap_range(1:end-1),CA1_mean_nday1y(1:end-1),CA1_nday1_sem(1:end-1),'or','linewidth',0.5,'MarkerFaceColor',usecolor2)
%scatter(CA1_nday1_x, CA1_nday1_y,'filled','MarkerFaceAlpha',0.01)


 [lap_range CA1_mean_nday2y CA1_nday2_sem]=cal_mean_sem(CA1_nday2_x, CA1_nday2_y);
 errorbar(lap_range(1:end-1),CA1_mean_nday2y(1:end-1),CA1_nday2_sem(1:end-1),'o','linewidth',0.5,'color', usecolor3,'MarkerFaceColor',usecolor3)
 %scatter(CA1_nday2_x, CA1_nday2_y,'filled','MarkerFaceAlpha',0.01)

legend({'f','nday1','nday2'})%,'nday2'})
ylim([-1.3 2])
title('CA1 f nday1 and nday2')
%% change to all point 
CA1_f_y(CA1_f_x==max(CA1_f_x))=[];
CA1_f_x(CA1_f_x==max(CA1_f_x))=[];
 CA1_f_model=fitlm(CA1_f_x, CA1_f_y);
plot(CA1_f_x,CA1_f_model.('Coefficients').('Estimate')(2).*CA1_f_x+ CA1_f_model.('Coefficients').('Estimate')(1),'color', usecolor1);
CA1_f_summary=anova(CA1_f_model,'summary');
CA1_f_p=CA1_f_summary.('pValue')(2);
CA1_R2=num2str(CA1_f_model.Rsquared.Ordinary)

%model0.('Coefficients').('Estimate')(2)
CA1_f_confidance=coefCI(CA1_f_model,0.001)
CA1_nday1_y(CA1_nday1_x==max(CA1_nday1_x))=[];
CA1_nday1_x(CA1_nday1_x==max(CA1_nday1_x))=[];
CA1_nday1_model=fitlm(CA1_nday1_x,CA1_nday1_y);
plot(CA1_nday1_x,CA1_nday1_model.('Coefficients').('Estimate')(2).*CA1_nday1_x+ CA1_nday1_model.('Coefficients').('Estimate')(1),'r');
CA1_nday1_confidance=coefCI(CA1_nday1_model,0.001)
CA1_nday1_summary=anova(CA1_nday1_model,'summary');

CA1_nday2_y(CA1_nday2_x==max(CA1_nday2_x))=[];
CA1_nday2_x(CA1_nday2_x==max(CA1_nday2_x))=[];
CA1_nday2_model=fitlm(CA1_nday2_x,CA1_nday2_y);
plot(CA1_nday2_x,CA1_nday2_model.('Coefficients').('Estimate')(2).*CA1_nday2_x+ CA1_nday2_model.('Coefficients').('Estimate')(1),'color', usecolor3);
CA1_nday2_confidance=coefCI(CA1_nday2_model,0.001)
CA1_nday2_summary=anova(CA1_nday2_model,'summary');
%%
subplot(1,2,2)
[lap_range CA3_mean_fy CA3_f_sem]=cal_mean_sem(CA3_f_x, CA3_f_y);
hold on;
%plot(range,mean_y,'o')
errorbar(lap_range(1:end-1),CA3_mean_fy(1:end-1),CA3_f_sem(1:end-1),'o','linewidth',0.5,'color', usecolor1,'MarkerFaceColor',usecolor1)

[lap_range CA3_mean_nday1y CA3_nday1_sem]=cal_mean_sem(CA3_nday1_x, CA3_nday1_y);
errorbar(lap_range(1:end-1),CA3_mean_nday1y(1:end-1),CA3_nday1_sem(1:end-1),'or','linewidth',0.5,'MarkerFaceColor',usecolor2)

% 
 [lap_range CA3_mean_nday2y CA3_nday2_sem]=cal_mean_sem(CA3_nday2_x, CA3_nday2_y);
 errorbar(lap_range(1:end-1),CA3_mean_nday2y(1:end-1),CA3_nday2_sem(1:end-1),'o','linewidth',0.5,'color', usecolor3,'MarkerFaceColor',usecolor3);
legend({'f','nday1','nday2'})%,'nday2'})
ylim([-1.3 2])
title('CA3 f nday1 and nday2')

CA3_f_y(CA3_f_x==max(CA3_f_x))=[];
CA3_f_x(CA3_f_x==max(CA3_f_x))=[];
CA3_f_model=fitlm(CA3_f_x, CA3_f_y);
plot(CA3_f_x,CA3_f_model.('Coefficients').('Estimate')(2).*CA3_f_x+ CA3_f_model.('Coefficients').('Estimate')(1),'color', usecolor1);
CA3_f_confidance=coefCI(CA3_f_model,0.001)
CA3_f_summary=anova(CA3_f_model,'summary');

CA3_nday1_y(CA3_nday1_x==max(CA3_nday1_x))=[];
CA3_nday1_x(CA3_nday1_x==max(CA3_nday1_x))=[];
CA3_nday1_model=fitlm(CA3_nday1_x, CA3_nday1_y);
plot(CA3_nday1_x,CA3_nday1_model.('Coefficients').('Estimate')(2).*CA3_nday1_x+ CA3_nday1_model.('Coefficients').('Estimate')(1),'r');
CA3_nday1_confidance=coefCI(CA3_nday1_model,0.001)
CA3_nday1_summary=anova(CA3_nday1_model,'summary');

CA3_nday2_y(CA3_nday2_x==max(CA3_nday2_x))=[];
CA3_nday2_x(CA3_nday2_x==max(CA3_nday2_x))=[];
CA3_nday2_model=fitlm(CA3_nday2_x, CA3_nday2_y);
plot(CA3_nday2_x,CA3_nday2_model.('Coefficients').('Estimate')(2).*CA3_nday2_x+ CA3_nday2_model.('Coefficients').('Estimate')(1),'color', usecolor3);
CA3_nday2_confidance=coefCI(CA3_nday2_model,0.001)
CA3_nday2_summary=anova(CA3_nday2_model,'summary');
set(gcf, 'Position',  [100, 100, 1200, 500])

%%
% figure;
% subplot(1,3,1);hold on;
% errorbar(lap_range(1:end-1),CA1_mean_fy(1:end-1),CA1_f_sem(1:end-1),'o','linewidth',0.5)
% errorbar(lap_range(1:end-1),CA3_mean_fy(1:end-1),CA3_f_sem(1:end-1),'o','linewidth',0.5)
% legend({'CA1 f','CA3 f'})
% ylim([-1.2 1.5])
% title('CA1 vs CA3 f ')
% 
% subplot(1,3,2);hold on;
% errorbar(lap_range(1:end-1),CA1_mean_nday1y(1:end-1),CA1_nday1_sem(1:end-1),'o','linewidth',0.5)
% errorbar(lap_range(1:end-1),CA3_mean_nday1y(1:end-1),CA3_nday1_sem(1:end-1),'o','linewidth',0.5)
% legend({'CA1 nday1','CA3 nday1'})
% ylim([-1.2 1.5])
% title('CA1 vs CA3 nday1 ')
% 
% subplot(1,3,3);hold on;
% errorbar(lap_range(1:end-1),CA1_mean_nday2y(1:end-1),CA1_nday2_sem(1:end-1),'o','linewidth',0.5)
% errorbar(lap_range(1:end-1),CA3_mean_nday2y(1:end-1),CA3_nday2_sem(1:end-1),'o','linewidth',0.5)
% legend({'CA1 nday2','CA3 nday2'})
% ylim([-1.2 1.5])
% title('CA1 vs CA3 nday2 ')
% set(gcf, 'Position',  [100, 100, 1600, 500])

%% 
figure
%model=fitlm(CA3_nday1_x, CA3_nday1_y);
errorbar(lap_range(1:end-1),CA1_mean_nday1y(1:end-1),CA1_nday1_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA1_nday1_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA1_nday1_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA1_nday1_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA1_nday1_confidance_2=coefCI(CA1_nday1_model,0.01)
%text(20,1,['CA1 R^2=' num2str(model1.('Rsquared').('Ordinary'))])

%CA1_nday2_model=fitlm(lap_range(1:end-1),CA3_mean_nday1y(1:end-1));
errorbar(lap_range(1:end-1),CA3_mean_nday1y(1:end-1),CA3_nday1_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA3_nday1_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA3_nday1_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA3_nday1_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA3_nday1_confidance_2=coefCI(CA3_nday1_model,0.01)
%text(20,0.5,['CA3 R^2=' num2str(model2.('Rsquared').('Ordinary'))])

title ('CA1 CA3 nday1, p<0.001')
legend({'CA1',' ','CA3',' ',})
box off
ylim([-1.5 2])
%ylim([-50 50]);

%%
figure
%CA1_nday1_model=fitlm(lap_range(1:end-1),CA1_mean_nday2y(1:end-1));
errorbar(lap_range(1:end-1),CA1_mean_nday2y(1:end-1),CA1_nday2_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA1_nday2_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA1_nday2_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA1_nday2_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA1_nday2_confidance_2=coefCI(CA1_nday2_model,0.01)
%text(20,1,['CA1 R^2=' num2str(model1.('Rsquared').('Ordinary'))])

%CA1_nday2_model=fitlm(lap_range(1:end-1),CA3_mean_nday2y(1:end-1));
errorbar(lap_range(1:end-1),CA3_mean_nday2y(1:end-1),CA3_nday2_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA3_nday2_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA3_nday2_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA3_nday2_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA3_nday2_confidance_2=coefCI(CA3_nday2_model,0.01)
%text(20,0.5,['CA3 R^2=' num2str(model2.('Rsquared').('Ordinary'))])

title ('CA1 CA3 nday2, p<0.001')
legend({'CA1',' ','CA3',' ',})
box off
ylim([-1.2 2])

%%
figure
%CA1_nday1_model=fitlm(lap_range(1:end-1),CA1_mean_fy(1:end-1));
errorbar(lap_range(1:end-1),CA1_mean_fy(1:end-1),CA1_f_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA1_f_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA1_f_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA1_f_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA1_f_confidance_2=coefCI(CA1_f_model,0.01)
%text(20,1,['CA1 R^2=' num2str(model1.('Rsquared').('Ordinary'))])
title ('CA1 f')
box off
ylim([-1.2 1.5])

figure;
%CA1_nday2_model=fitlm(lap_range(1:end-1),CA3_mean_fy(1:end-1));
errorbar(lap_range(1:end-1),CA3_mean_fy(1:end-1),CA3_f_sem(1:end-1),'o','linewidth',0.5)
hold on;
a=CA3_f_model.('Coefficients').('Estimate')(2);
plot(lap_range(1:end-1),CA3_f_model.('Coefficients').('Estimate')(2).*lap_range(1:end-1)+ CA3_f_model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
CA3_f_confidance_2=coefCI(CA3_f_model,0.01)
%text(20,0.5,['CA3 R^2=' num2str(model2.('Rsquared').('Ordinary'))])

title ('ca3 f')
%legend({'CA1',' ','CA3',' ',})
box off
ylim([-1.2 1.5])