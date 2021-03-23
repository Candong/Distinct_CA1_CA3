clear all; close all;
% for now the slope is calculated by onstlap and later com differences, it
% could also include the earlier activities before 
% 
[pf_files, temp]=uigetfile('*.mat', 'Chose PF files to load:','MultiSelect','on');

f_count=1;
n_count=1;
for f=1:size(pf_files,2)
    if contains(pf_files{f},'_f_')
       f_pf_files{f_count}= pf_files{f};
       f_count=f_count+1;
    elseif contains(pf_files{f},'_n_')
        n_pf_files{n_count}= pf_files{f};
        n_count=n_count+1;
    
    end
end
f_pf_files=sort(f_pf_files);
n_pf_files=sort(n_pf_files);


window=1;
onsetlap=1;
step=5;
fit_lapnum=[];
binsize=6;
%%
f_slope=[];
f_all_start_lap=[];
f_COM_start=[];
f_COM_end=[];
f_COM_alllaps=[];
f_total_meantrans=[];
f_onset_deltaCOM=[];
f_Rsquare=[];
%f_all_deltaCOM=[];
f_pvalue=[];
for f=1:size(f_pf_files,2)
load(f_pf_files{f});
[slope pvalue all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM all_deltaCOM pf_id Rsquare]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
f_slope=[f_slope slope];
f_all_start_lap=[f_all_start_lap all_start_lap];
f_COM_start=[f_COM_start COM_start];
f_COM_end=[f_COM_end COM_end];
f_COM_alllaps=[f_COM_alllaps COM_alllaps];
f_pf_id{f}=pf_id;
f_mean_trans{f}=mean_trans(pf_id,:);
f_total_meantrans = [f_total_meantrans; mean_trans(pf_id,:)]; 
f_all_deltaCOM {f}=all_deltaCOM;
f_startlap{f}=all_start_lap;
f_onset_deltaCOM=[f_onset_deltaCOM onset_deltaCOM];
f_Rsquare=[f_Rsquare Rsquare];
f_pvalue=[f_pvalue pvalue];
end
%%
yrange=[-0.8 0.8];
plot_real_abs_fit(f_COM_alllaps*binsize,f_slope*binsize,f_pvalue,yrange*binsize,10*binsize)
legend({'sig forward','sig backward'})
title('CA3 f')

%%
yrange=[-0.8 0.8];
plot_real_abs_fit(f_all_start_lap(f_all_start_lap<=30),f_slope(f_all_start_lap<=30)*binsize,f_pvalue(f_all_start_lap<=30),yrange*binsize,10)
legend({'sig forward','sig backward'})
title('CA3 f')
xlabel('PF onset lap')
xlim([0 30])
%%
%[pf_files, temp]=uigetfile('*.mat', 'Chose PF files to load:','MultiSelect','on');
n_slope=[];
n_all_start_lap=[];
n_COM_start=[];
n_COM_end=[];
n_COM_alllaps=[];
n_total_meantrans=[];
n_all_deltaCOM=[];
n_onset_deltaCOM=[];
n_Rsquare=[];
n_pvalue=[];
for f=1:size(n_pf_files,2)
load(n_pf_files{f});
[slope pvalue all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM all_deltaCOM pf_id Rsquare]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
n_slope=[n_slope slope];
n_all_start_lap=[n_all_start_lap all_start_lap];
n_COM_start=[n_COM_start COM_start];
n_COM_end=[n_COM_end COM_end];
n_COM_alllaps=[n_COM_alllaps COM_alllaps];
f_n_pf_id=f_pf_id{f};
n_mean_trans{f}=mean_trans(f_n_pf_id,:);
n_total_meantrans = [n_total_meantrans; mean_trans(f_n_pf_id,:)]; 
n_all_deltaCOM{f}= all_deltaCOM;
n_startlap{f}=all_start_lap;
n_onset_deltaCOM=[n_onset_deltaCOM onset_deltaCOM];
n_Rsquare =[n_Rsquare Rsquare];
n_pvalue = [n_pvalue pvalue];
end
%%
% figure;
% plot(n_slope*binsize,n_Rsquare,'o')
% xlabel('slope')
% ylabel('Rsquare')
% title('novel')
% box off

%%
%significant
figure;
hold on
plot(n_slope(n_pvalue<0.05 &  n_slope<0)*binsize,n_Rsquare(n_pvalue<0.05 &  n_slope<0),'o')
plot(n_slope(n_pvalue<0.05 &  n_slope>0)*binsize,n_Rsquare(n_pvalue<0.05&  n_slope>0),'o')

plot(n_slope(n_pvalue>=0.05)*binsize,n_Rsquare(n_pvalue>=0.05),'o','Color',[0,0,0]+0.6)
legend({'significant','insignificant'})
xlabel('slope')
ylabel('Rsquare')
title('novel')
box off
%%
yrange=[-0.8 0.8];

plot_real_abs_fit(n_COM_alllaps*binsize,n_slope*binsize,n_pvalue,yrange*binsize,10*binsize)
legend({'sig forward','sig backward'})
title('CA1 nday2')
%%
yrange=[-0.8 0.8];
plot_real_abs_fit(n_all_start_lap(n_all_start_lap<=30),n_slope(n_all_start_lap<=30)*binsize,n_pvalue(n_all_start_lap<=30),yrange*binsize,10)
legend({'sig forward','sig backward'})
title('CA1 nady2')
xlabel('PF onset lap')
xlim([0 30])
%%
% sig percentage by lap
sig_per_bylap=[];
lap=1:15;
for i=lap
    cur_p=n_pvalue(n_all_start_lap==i);
    cur_percentage=sum(cur_p<0.05)/length(cur_p);
    sig_per_bylap=[sig_per_bylap cur_percentage];
    
    
    
end
figure
bar(sig_per_bylap)


figure;
hold on;
plot(sig_per_bylap,'o')
model = fitlm(lap,sig_per_bylap);
plot(lap,model.('Coefficients').('Estimate')(2).*lap+ model.('Coefficients').('Estimate')(1));
model_summary=anova(model,'summary')
text(1, 0.5, ['p =' num2str(model_summary.('pValue')(2))])
%%
figure;
slope_propotion=[sum(n_pvalue<0.05  & n_slope<0) sum(n_pvalue<0.05  & n_slope>0) sum(n_pvalue>0.05)];
labels={'sig back','sig fore','non sig'};
bar(1,slope_propotion./length(n_pvalue))
legend(labels)
title('CA1  n')
%%
%%f_slope
 fn_cor=diag(corr(f_total_meantrans', n_total_meantrans'));
 yrange=[-0.8 0.8];

plot_real_abs_fit(fn_cor', f_slope*binsize,f_pvalue,yrange*binsize,-0.2)

% figure;hold on;
% plot(fn_cor,f_slope,'o');
% model_f=fitlm(fn_cor, f_slope);
% plot(fn_cor,model_f.('Coefficients').('Estimate')(2).*fn_cor+ model_f.('Coefficients').('Estimate')(1));
% f_CI=coefCI(model_f,0.05)
% title(['slope = ' num2str(model_f.('Coefficients').('Estimate')(2))*6])
%subplot(1,2,1)
xlabel('novel day1 day2 correlation')
ylabel('slope')
legend({'sig forward','sig backward'})
title('CA3 slope vs twodaycor')

% subplot(1,2,2)
% xlabel('novel day1 day2 correlation')
% ylabel('slope')
% title('abs slope vs twodaycor')
% text(0.5,0.1,['intercept ' num2str(f_CI(1,:))])
% text(0.5,0.12,['slope ' num2str(f_CI(2,:))])
%%
high_cor_slope=f_slope(fn_cor>=0.5);
low_cor_slope=f_slope(fn_cor<0.5);

[h p]=ttest2(high_cor_slope, low_cor_slope)


%% by onset lap delta COM
% keep_id=~isnan(f_onset_deltaCOM);
% new_f_all_start_lap=f_all_start_lap(keep_id);
% new_f_onset_deltaCOM=f_onset_deltaCOM(keep_id);
% figure; hold on;
% %plot(new_f_all_start_lap,new_f_onset_deltaCOM,'o')
% title('familiar')
% xlabel('onset lap')
% ylabel('by lap delta COM')
% [f_lap_range f_mean_onset_deltaCOM f_sem_onset_deltaCOM]=cal_mean_sem(new_f_all_start_lap,new_f_onset_deltaCOM);
% %figure
% errorbar(f_lap_range,f_mean_onset_deltaCOM,f_sem_onset_deltaCOM,'or','linewidth',1.5)
% xlim([0 10])

keep_id=~isnan(n_onset_deltaCOM);
new_n_all_start_lap=n_all_start_lap(keep_id);
new_n_onset_deltaCOM=n_onset_deltaCOM(keep_id);
% figure; hold on;
% %plot(new_n_all_start_lap,new_n_onset_deltaCOM,'o')
% title('novel')
% xlabel('onset lap')
% ylabel('by lap delta COM')
% [n_lap_range n_mean_onset_deltaCOM n_sem_onset_deltaCOM]=cal_mean_sem(new_n_all_start_lap,new_n_onset_deltaCOM);
% %figure
% errorbar(n_lap_range,n_mean_onset_deltaCOM,n_sem_onset_deltaCOM,'or','linewidth',1.5)
% 
% xlim([0 10])
%% plot boxplot version of onsetlap deltaCOM vs onsetlap
figure;
G=[];
COM_v=[];
order={};
for i=1:10
    cur_id = new_n_all_start_lap==i;
    cur_onset_deltaCOM=new_n_onset_deltaCOM(cur_id);
    for g=1:sum(cur_id)
        if i<10
           
            G=[G; [num2str(i) ' ']];
        else 
            G=[G; num2str(i)];
        end
    end
    COM_v=[COM_v cur_onset_deltaCOM];
    
    order{i,1}=char(num2str(i));
%     figure; histogram(cur_onset_deltaCOM,-49.5:1:49.5)
end
G=cellstr(G);
figure
violinplot(COM_v*6,G,'GroupOrder',order);
ylabel('delta COM nextt lap-onset  lap (cm)')
xlabel('onsetlap')
title('CA1 deltaCOM nday2')

%%
figure; 
%subplot(1,2,1)
hold on
histogram(new_n_onset_deltaCOM(new_n_all_start_lap==1)*6,(-49.5:1:49.5)*6,'Normalization','Probability');
histogram(new_n_onset_deltaCOM(new_n_all_start_lap>1)*6,(-49.5:1:49.5)*6,'Normalization','Probability');
xlim([-180 180])
legend({'instant PF','delayed  PF'})
title('delta COM onsetlap and later lap')
% subplot(1,2,2)
% hold on;
% cdfplot(new_n_onset_deltaCOM(new_n_all_start_lap==1))
% cdfplot(new_n_onset_deltaCOM(new_n_all_start_lap>1))
% legend({'instant PF','delayed  PF'})
% title('cdf od instant and delayedPF onset lap deltaCOM')
% grid off


%set(gcf, 'Position',  [100, 100, 1200, 500])
%p_deltacom=ranksum(new_n_onset_deltaCOM(new_n_all_start_lap==1),new_n_onset_deltaCOM(new_n_all_start_lap>1))

[h p_deltacom]=kstest2(new_n_onset_deltaCOM(new_n_all_start_lap==1),new_n_onset_deltaCOM(new_n_all_start_lap>1))
text(80,0.2,['kstest2 p= ' num2str(p_deltacom)])
%% by lap delta COM
% [f_delta_com f_lap]= shapecell2vector_COM_lap(f_all_deltaCOM);
% [n_delta_com n_lap]= shapecell2vector_COM_lap(n_all_deltaCOM);
% figure; hold on;
% %plot(f_lap,f_delta_com,'o')
% title('familiar')
% xlabel('lap')
% ylabel('by lap delta COM')
% [f_lap_range f_mean_dcom f_sem_dcom]=cal_mean_sem(f_lap,f_delta_com);
% errorbar(f_lap_range,f_mean_dcom,f_sem_dcom,'or','linewidth',1.5)
% 
% figure; hold on;
% %plot(n_lap,n_delta_com,'o')
% title('novel')
% xlabel('lap')
% ylabel('by lap delta COM')
% [n_lap_range n_mean_dcom n_sem_dcom]=cal_mean_sem(n_lap,n_delta_com);
% errorbar(n_lap_range,n_mean_dcom,n_sem_dcom,'or','linewidth',1.5)
%%
% function []=plot_real_abs_fit(n_COM_alllaps,n_slope,yrange)
% 
% figure;
% subplot(1,2,1)
% hold on
% plot(n_COM_alllaps,n_slope,'o')
% title('novel')
% xlabel('PF position')
% ylabel('slope')
% model_n=fitlm(n_COM_alllaps,n_slope);
% plot(n_COM_alllaps,model_n.('Coefficients').('Estimate')(2).*n_COM_alllaps+ model_n.('Coefficients').('Estimate')(1));
% %title(['slope = ' num2str(model_n.('Coefficients').('Estimate')(2)*6)])
% text(25,yrange12)+0.3,['intercept ' num2str(model_n.('Coefficients').('Estimate')(1)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(1)) ])
% text(25,yrange(1)+0.2,['slope ' num2str(model_n.('Coefficients').('Estimate')(2)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(2)) ])
% text(25,yrange(1)+0.1,['R-square ' num2str(model_n.Rsquared.Ordinary)]);
% title ('F real value slope vs position')
% ylim(yrange)
% 
% subplot(1,2,2)
% hold on;
% plot(n_COM_alllaps,abs(n_slope),'o')
% title('novel')
% xlabel('PF position')
% ylabel('slope')
% model_n=fitlm(n_COM_alllaps,abs(n_slope));
% plot(n_COM_alllaps,model_n.('Coefficients').('Estimate')(2).*n_COM_alllaps+ model_n.('Coefficients').('Estimate')(1));
% %title(['slope = ' num2str(model_f.('Coefficients').('Estimate')(2)*6)])
% text(25,yrange(1)+0.3,['intercept ' num2str(model_n.('Coefficients').('Estimate')(1)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(1)) ])
% text(25,yrange(1)+0.2,['slope ' num2str(model_n.('Coefficients').('Estimate')(2)) 'p-value ' num2str(model_n.('Coefficients').('pValue')(2)) ])
% text(25,yrange(1)+0.1,['R-square ' num2str(model_n.Rsquared.Ordinary)]);title('F abs value slope vs posotopn')
% ylim(yrange)
% set(gcf, 'Position',  [100, 100, 1200, 500])
% 
% 
% 
% 
% 
% end



%%
% function [f_delta_com f_lap]= shapecell2vector_COM_lap(f_all_deltaCOM)
% f_delta_com=[];
% f_lap=[];
% for f=1:size(f_all_deltaCOM,2)
%     cur_f_dcom=f_all_deltaCOM{f};
%     cur_f_lap=ones(size(cur_f_dcom,1),size(cur_f_dcom,2)).*(1:size(cur_f_dcom,2));
%     cur_f_dcom=reshape(cur_f_dcom,1,[]);
%     cur_f_lap=reshape(cur_f_lap,1,[]);
%     f_remove_id=isnan(cur_f_dcom);
%     cur_f_dcom(f_remove_id)=[];
%     cur_f_lap(f_remove_id)=[];
%     f_delta_com=[f_delta_com cur_f_dcom];
%     f_lap=[f_lap cur_f_lap];
%     %cur_n_dcom=n_all_deltaCOM{f};
%     
%     
%     
% end
% end


%%
% function [slope all_start_lap COM_start COM_end COM_alllaps]=caculate_shift_parameters(sig_PFs,window,onsetlap,step)
% delta_COM=[];
% Rsquare=[];
% slope=[];
% Rsquare2=[];
% slope=[];
% var_y=[];
% var_y2=[];
% 
% [COM all_start_lap]=calculate_and_linreg_eachCOM_by_window(sig_PFs,window,onsetlap);
% [COM_start COM_end COM_alllaps]=calculate_start_end_alllaps_COM(sig_PFs,step,onsetlap);
% 
% % total_mean_var=[];
% %figure; hold on;
%     for i=1:size(COM,1)
%         cur_deltaCOM= COM(i,:)-COM(i,all_start_lap(i)); 
%         delta_COM=[delta_COM; cur_deltaCOM ];
%         % figure;hold on;
%         % plot(1:size(COM,2), cur_deltaCOM,'o');
%         % ylim([-10 10])
% 
%         x=1:size(COM,2);
%         x(isnan(cur_deltaCOM))=[];
%         y=cur_deltaCOM(~isnan(cur_deltaCOM));
%         if onsetlap==0 %start regression from first lap with activity
% %             cur_var_y=var(y);
% %             var_y=[var_y cur_var_y];
%             model=fitlm(x, y);
%             % plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
%             % xlim([1 size(COM,2)]);
% 
%             slope=[slope model.('Coefficients').('Estimate')(2)];
%             cur_r_squared=model.Rsquared.Ordinary;
%             Rsquare=[Rsquare cur_r_squared];
% 
%         else
%             x2=x(x>=all_start_lap(i));
%             y2=y(x>=all_start_lap(i));
% %             cur_var_y2=var(y2);
% %             var_y2=[var_y2 cur_var_y2];
%             model=fitlm(x2,y2);
%             slope=[slope model.('Coefficients').('Estimate')(2)];
%             cur_r_squared2=model.Rsquared.Ordinary;
%             Rsquare2=[Rsquare2 model.Rsquared.Ordinary];
%             mean_var=windowed_var(y2,3);
%             total_mean_var=[total_mean_var mean_var];
% 
%         end
% 
% 
% 
% 
% 
%     end
% end



