clear all; %close all;
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
fit_lapnum=[1 5];
%%
%fit to lap 1:5
f_slope=[];
f_all_start_lap=[];
f_COM_start=[];
f_COM_end=[];
f_COM_alllaps=[];

for f=1:size(f_pf_files,2)
load(f_pf_files{f});
[slope all_start_lap COM_start COM_end COM_alllaps]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
f_slope=[f_slope slope];
f_all_start_lap=[f_all_start_lap all_start_lap];
f_COM_start=[f_COM_start COM_start];
f_COM_end=[f_COM_end COM_end];
f_COM_alllaps=[f_COM_alllaps COM_alllaps];




end

% figure;plot(f_COM_start,f_slope,'o')
% figure;plot(f_COM_end,f_slope,'o')
remove_id=f_slope==-100;
f_slope_copy=f_slope;
f_all_start_lap_copy= f_all_start_lap;
f_COM_alllaps(remove_id)=[];
f_slope(remove_id)=[];

% figure;plot(f_COM_alllaps,f_slope,'o')
% title('familiar')
% xlabel('PF position')
% ylabel('slope')
% 
% 
% f_all_start_lap(remove_id)=[];
% 
% figure; hold on;
% plot(f_all_start_lap,f_slope,'o')
% title('familiar')
% xlabel('start lap')
% ylabel('slope')
% [lap_range f_mean_slope f_sem_slope]=cal_mean_sem(f_all_start_lap, f_slope);
% errorbar(lap_range,f_mean_slope,f_sem_slope,'or','linewidth',1)

%%
%fit to lap 6-10
f_slope2=[];
fit_lapnum=[6 10];
f_all_start_lap2=[];
f_COM_start2=[];
f_COM_end2=[];
f_COM_alllaps2=[];
f_onset_deltaCOM=[];
for f=1:size(f_pf_files,2)
load(f_pf_files{f});
[slope all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
f_slope2=[f_slope2 slope];
f_all_start_lap2=[f_all_start_lap2 all_start_lap];
f_COM_start2=[f_COM_start2 COM_start];
f_COM_end2=[f_COM_end2 COM_end];
f_COM_alllaps2=[f_COM_alllaps2 COM_alllaps];
f_onset_deltaCOM=[f_onset_deltaCOM onset_deltaCOM];



end

% figure;plot(f_COM_start,f_slope,'o')
% figure;plot(f_COM_end,f_slope,'o')
remove_id=f_slope2==-100;
f_slope2_copy=f_slope2;
f_all_start_lap2_copy=f_all_start_lap2;
f_COM_alllaps2(remove_id)=[];
f_slope2(remove_id)=[];

% figure;plot(f_COM_alllaps2,f_slope2,'o')
% title('familiar')
% xlabel('PF position')
% ylabel('slope')
% 
% 
% f_all_start_lap2(remove_id)=[];
% 
% figure; hold on;
% plot(f_all_start_lap2,f_slope2,'o')
% title('familiar')
% xlabel('start lap')
% ylabel('slope')
% [lap_range f_mean_slope f_sem_slope]=cal_mean_sem(f_all_start_lap2, f_slope2);
% errorbar(lap_range,f_mean_slope,f_sem_slope,'or','linewidth',1)



%%
%[pf_files, temp]=uigetfile('*.mat', 'Chose PF files to load:','MultiSelect','on');
n_slope=[];
n_all_start_lap=[];
n_COM_start=[];
n_COM_end=[];
n_COM_alllaps=[];
fit_lapnum=[1 5];
n_onset_deltaCOM=[];
for f=1:size(n_pf_files,2)
load(n_pf_files{f});
[slope all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
n_slope=[n_slope slope];
n_all_start_lap=[n_all_start_lap all_start_lap];
n_COM_start=[n_COM_start COM_start];
n_COM_end=[n_COM_end COM_end];
n_COM_alllaps=[n_COM_alllaps COM_alllaps];
n_onset_deltaCOM=[n_onset_deltaCOM onset_deltaCOM];





end

remove_id=n_slope==-100;
n_slope_copy=n_slope;
n_all_start_lap_copy=n_all_start_lap;
n_COM_alllaps(remove_id)=[];
n_slope(remove_id)=[];
% figure;plot(n_COM_start,n_slope,'o')
% figure;plot(n_COM_end,n_slope,'o')
% figure;plot(n_COM_alllaps,n_slope,'o')
% title('novel')
% xlabel('PF position')
% ylabel('slope')
% 
% n_all_start_lap(remove_id)=[];
% 
% figure; hold on;
% plot(n_all_start_lap,n_slope,'o')
% title('novel')
% xlabel('start lap')
% ylabel('slope')
% [n_lap_range n_mean_slope n_sem_slope]=cal_mean_sem(n_all_start_lap, n_slope);
% errorbar(n_lap_range,n_mean_slope,n_sem_slope,'or','linewidth',1)

%%
n_slope2=[];
n_all_start_lap2=[];
n_COM_start2=[];
n_COM_end2=[];
n_COM_alllaps2=[];
fit_lapnum=[6 10];

for f=1:size(n_pf_files,2)
load(n_pf_files{f});
[slope all_start_lap COM_start COM_end COM_alllaps]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
n_slope2=[n_slope2 slope];
n_all_start_lap2=[n_all_start_lap2 all_start_lap];
n_COM_start2=[n_COM_start2 COM_start];
n_COM_end2=[n_COM_end2 COM_end];
n_COM_alllaps2=[n_COM_alllaps2 COM_alllaps];





end

remove_id=n_slope2==-100;
n_slope2_copy=n_slope2;
n_all_start_lap2_copy=n_all_start_lap2;
n_COM_alllaps2(remove_id)=[];
n_slope2(remove_id)=[];
% figure;plot(n_COM_start,n_slope,'o')
% figure;plot(n_COM_end,n_slope,'o')
figure;plot(n_COM_alllaps2,n_slope2,'o')
title('novel')
xlabel('PF position')
ylabel('slope')

n_all_start_lap2(remove_id)=[];

% figure; hold on;
% plot(n_all_start_lap2,n_slope2,'o')
% title('novel')
% xlabel('start lap')
% ylabel('slope')
% [n_lap_range n_mean_slope n_sem_slope]=cal_mean_sem(n_all_start_lap2, n_slope2);
% errorbar(n_lap_range,n_mean_slope,n_sem_slope,'or','linewidth',1)
%%
f_plot_id=f_slope_copy~=-100 & f_slope2_copy~=-100 & f_all_start_lap_copy==1;
figure; hold on;
plot(ones(1,sum(f_plot_id)),f_slope_copy(f_plot_id),'o');
plot(ones(1,sum(f_plot_id))*2,f_slope2_copy(f_plot_id),'o');
plot([ones(1,sum(f_plot_id)); ones(1,sum(f_plot_id))*2],[f_slope_copy(f_plot_id);f_slope2_copy(f_plot_id)])
xlim([0 3])
title('familiar')
[f_h f_p]=ttest(f_slope_copy(f_plot_id),f_slope2_copy(f_plot_id))

figure;plot(f_slope_copy(f_plot_id),f_slope2_copy(f_plot_id),'o')
hold on;plot([-5 4],[-5 4])
xlabel('1 to 5 lap fit')
ylabel('6 to 10 lap fit')

n_plot_id=n_slope_copy~=-100 & n_slope2_copy~=-100 & n_all_start_lap_copy==1;
figure; hold on;
plot(ones(1,sum(n_plot_id)),n_slope_copy(n_plot_id),'o');
plot(ones(1,sum(n_plot_id))*2,n_slope2_copy(n_plot_id),'o');
plot([ones(1,sum(n_plot_id)); ones(1,sum(n_plot_id))*2],[n_slope_copy(n_plot_id);n_slope2_copy(n_plot_id)])
xlim([0 3])
title('novel')
[n_h n_p]=ttest(n_slope_copy(n_plot_id),n_slope2_copy(n_plot_id))
figure;plot(n_slope_copy(n_plot_id),n_slope2_copy(n_plot_id),'o')
hold on;plot([-5 4],[-5 4])
xlabel('1 to 5 lap fit')
ylabel('6 to 10 lap fit')

%%
figure;histogram(f_onset_deltaCOM(f_all_start_lap_copy==1),-24.5:1:24.5,'Normalization','probability')
hold on; histogram(f_onset_deltaCOM(f_all_start_lap_copy>1),-24.5:1:24.5,'Normalization','probability')
title('familiar onsetlap+1 - onsetlap COM')
legend({'instant','delayed'})



figure;histogram(n_onset_deltaCOM(n_all_start_lap_copy==1),-24.5:1:24.5,'Normalization','probability')
hold on; histogram(n_onset_deltaCOM(n_all_start_lap_copy>1),-24.5:1:24.5,'Normalization','probability')
title('novel onsetlap+1 - onsetlap COM')
legend({'instant','delayed'})
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



