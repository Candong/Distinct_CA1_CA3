clear all; close all;

%%
[v_info, temp]=uigetfile('*.mat', 'Chose v_info files to load:','MultiSelect','on');



total_f_s=[];
total_f_v=[];
total_mean_f_s=[];
total_mean_f_v=[];
total_n_s=[];
total_n_v=[];
total_mean_n_s=[];
total_mean_n_v=[];
total_nday2_s=[];
total_nday2_v=[];
total_mean_nday2_s=[];
total_mean_nday2_v=[];

for f=1:size(v_info,2)
    load(v_info{f});
    day1=~isempty(fieldnames(f1_v));
    day2=~isempty(fieldnames(f2_v));
    nday2=~isempty(fieldnames(nday2_v));
    
    day1_id(f)=day1;
    day2_id(f)=day2;
    nday2_id(f)=nday2;
    
    if day1
        [cur_f_v cur_f_s mean_cur_f_v mean_cur_f_s]=align_shift_v(f1_v,day1_f_COM_diff);
        [cur_n_v cur_n_s mean_cur_n_v mean_cur_n_s]=align_shift_v(n1_v,day1_n_COM_diff);
        total_f_s=[total_f_s cur_f_s];
        total_f_v=[total_f_v cur_f_v];
        total_mean_f_s=[total_mean_f_s total_mean_f_s];
        total_mean_f_v=[total_mean_f_v total_mean_f_v];
        total_n_s=[total_n_s cur_n_s];
        total_n_v=[total_n_v cur_n_v];
        total_mean_n_s=[total_mean_n_s total_mean_n_s];
        total_mean_n_v=[total_mean_n_v total_mean_n_v];        
    end
    
    if day2
        [cur_f_v cur_f_s mean_cur_f_v mean_cur_f_s]=align_shift_v(f2_v,day2_f_COM_diff);
        [cur_n_v cur_n_s mean_cur_n_v mean_cur_n_s]=align_shift_v(n2_v,day2_n_COM_diff);
        total_f_s=[total_f_s cur_f_s];
        total_f_v=[total_f_v cur_f_v];
        total_mean_f_s=[total_mean_f_s total_mean_f_s];
        total_mean_f_v=[total_mean_f_v total_mean_f_v];
        total_n_s=[total_n_s cur_n_s];
        total_n_v=[total_n_v cur_n_v];
        total_mean_n_s=[total_mean_n_s total_mean_n_s];
        total_mean_n_v=[total_mean_n_v total_mean_n_v];         
    end
    
    if nday2
        [cur_nday2_v cur_nday2_s mean_cur_nday2_v mean_cur_nday2_s]=align_shift_v(nday2_v,nday2_n_COM_diff);
        total_nday2_s=[total_nday2_s cur_nday2_s];
        total_nday2_v=[total_nday2_v cur_nday2_v];
        total_mean_nday2_s=[total_mean_nday2_s total_mean_nday2_s];
        total_mean_nday2_v=[total_mean_nday2_v total_mean_nday2_v];          
    end
    
    
    
    
end
%%
figure; hold on;
plot(total_f_v,total_f_s,'o')
title('familiar')
ylabel('lap COM difference/ bin')
xlabel('velocity/cm')
model_f=fitlm(total_f_v, total_f_s);
plot(total_f_v,model_f.('Coefficients').('Estimate')(2).*total_f_v+ model_f.('Coefficients').('Estimate')(1));
f_CI=coefCI(model_f,0.05);
f_summary=anova(model_f,'summary');
text(25,-10,['intercept ' num2str(f_CI(1,:))])
text(25,-12,['slope ' num2str(f_CI(2,:))])

figure; hold on;
plot(total_n_v,total_n_s,'o')
title('novel')
ylabel('lap COM difference/ bin')
xlabel('velocity/cm')
model_n=fitlm(total_n_v, total_n_s);
plot(total_n_v,model_n.('Coefficients').('Estimate')(2).*total_n_v+ model_n.('Coefficients').('Estimate')(1));
n_CI=coefCI(model_n,0.05);
n_summary=anova(model_n,'summary');
text(25,-10,['intercept ' num2str(n_CI(1,:))])
text(25,-12,['slope ' num2str(n_CI(2,:))])

figure; hold on;
plot(total_nday2_v,total_nday2_s,'o')
title('nday2')
ylabel('lap COM difference/ bin')
xlabel('velocity/cm')
model_nday2=fitlm(total_nday2_v, total_nday2_s);
plot(total_nday2_v,model_nday2.('Coefficients').('Estimate')(2).*total_nday2_v+ model_n.('Coefficients').('Estimate')(1));
nday2_CI=coefCI(model_nday2,0.05);
nday2_summary=anova(model_nday2,'summary');
text(25,-10,['intercept ' num2str(nday2_CI(1,:))])
text(25,-12,['slope ' num2str(nday2_CI(2,:))])