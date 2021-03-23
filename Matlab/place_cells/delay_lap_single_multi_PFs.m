clear all; close all;
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;

if ~isa(behavior,'cell') 
  behavior={behavior};  
    
end 
for f=1:size(behavior,2)
    if contains(behavior{f},'_f_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);



PF_center=[];
PF_mean=[];
PF_count=1;
window=6;
threshold=3;
PF_start_lap=[];
range=1:30;
f_total_act_lap=[];
f_total_firstlap_act=[];
f_single_start_lap=[];
f_multi_start_lap=[];
f_count_single=0;
f_count_total=0;
f_count_multi=0;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
single_id= find(number_of_PFs==1);
multi_id=find(number_of_PFs >1);
f_count_total=f_count_total+length(single_id)+length(multi_id);
f_count_single=f_count_single+length(single_id);
f_count_multi=f_count_multi+length(multi_id);
for i=1:size(sig_PFs,2) 
    %f_count_total= f_count_total+1;
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs_with_noise{j,i};
            binmean_PF_temp=sig_PFs{j,i};            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2); 
            start_bin=PF_start_bins(j,i);
            end_bin=PF_end_bins(j,i);
            max_lap=30;%size(binmean_PF_temp,2);
%             PF_center(PF_count)=find(mean_PF==max(mean_PF));
            start_lap=0;lap=1;
            start_lap=find_delaylap(binmean_PF_temp,start_bin,end_bin,window,threshold,max_lap);
            if start_lap>1
               act_lap=sum(mean( binmean_PF_temp(start_bin:end_bin,:),2)>0);
               firstlap_act=sum( binmean_PF_temp(start_bin:end_bin,1));
               f_total_act_lap=[f_total_act_lap act_lap];
               f_total_firstlap_act=[f_total_firstlap_act firstlap_act];
               if firstlap_act>0
                   %figure;imagesc(binmean_PF_temp');
               end  
            end

            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
            if ismember(i, single_id)
                f_single_start_lap=[f_single_start_lap start_lap];
                
            elseif ismember(i, multi_id)
                f_multi_start_lap=[f_multi_start_lap start_lap];
                
            end
                
        end
    
    end
                   
end
f_PF_id=find(isnan(number_of_PFs));



end
f_PF_start=PF_start_lap;
f_distribution=histc(PF_start_lap,range);

%ylim([0 350]);

 display(['n day1 percent' ]);
 display(sum(f_total_firstlap_act>0)/length(f_total_firstlap_act));
%%
figure;
subplot(1,2,1)
hold on;
histogram(f_single_start_lap,1:30,'Normalization','Probability')
histogram(f_multi_start_lap,1:30,'Normalization','Probability')
title('Pf onset lap familiar')
legend({'single','multi'})
text(10,0.2, [ 'signle PF cell ratio ' num2str(f_count_single/f_count_total)])
text(10,0.3, [ 'multi PF cell ratio ' num2str(f_count_multi/f_count_total)])


subplot(1,2,2)
hold on;
cdfplot(f_single_start_lap)
cdfplot(f_multi_start_lap)
title('Pf onset lap')
legend({'single','multi'})
[h p]=kstest2(f_single_start_lap,f_multi_start_lap);
text(10,0.2, [ 'kstest2 p=' num2str(p)])




%% 
file_count=1;

if ~isa(behavior,'cell') 
  behavior={behavior};  
    
end 
for f=1:size(behavior,2)
    if contains(behavior{f},'_n_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end


behavior_filepaths=sort(behavior_filepaths);PF_center=[];


PF_mean=[];
PF_count=1;
window=6;
threshold=3;
PF_start_lap=[];
range=1:30;
n_total_act_lap=[];
n_total_firstlap_act=[];
n_single_start_lap=[];
n_multi_start_lap=[];
n_count_single=0;
n_count_total=0;
n_count_multi=0;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
single_id= find(number_of_PFs==1);
multi_id=find(number_of_PFs >1);
n_count_total=n_count_total+length(single_id)+length(multi_id);
n_count_single=n_count_single+length(single_id);
n_count_multi=n_count_multi+length(multi_id);
for i=1:size(sig_PFs,2)
    for j=1%:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs_with_noise{j,i};
            binmean_PF_temp=sig_PFs{j,i};            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2);            
            start_bin=PF_start_bins(j,i);
            end_bin=PF_end_bins(j,i);
            max_lap=30;%size(binmean_PF_temp,2);
%             PF_center(PF_count)=find(mean_PF==max(mean_PF));
            start_lap=0;lap=1;
            start_lap=find_delaylap(binmean_PF_temp,start_bin,end_bin,window,threshold,max_lap);
            if start_lap>1
               act_lap=sum(mean( binmean_PF_temp(start_bin:end_bin,:),2)>0);
               firstlap_act=sum( binmean_PF_temp(start_bin:end_bin,1));
%                if firstlap_act>0
%                    figure;imagesc(binmean_PF_temp');
%                end
               n_total_act_lap=[n_total_act_lap act_lap];
               n_total_firstlap_act=[n_total_firstlap_act firstlap_act];
              
            end

            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
            if ismember(i, single_id)
                n_single_start_lap=[n_single_start_lap start_lap];
            elseif ismember(i, multi_id)
                n_multi_start_lap=[n_multi_start_lap start_lap];
                
            end
        end
    
    end
                   
end


end
 display(['n day2 percent' ]);
 display(sum(n_total_firstlap_act>0)/length(n_total_firstlap_act));
 n_PF_start=PF_start_lap;
 
 %%
 figure;
subplot(1,2,1)
hold on;
histogram(n_single_start_lap,1:30,'Normalization','Probability')
histogram(n_multi_start_lap,1:30,'Normalization','Probability')
title('Pf onset lap novel')
legend({'single','multi'})
text(10,0.05, [ 'signle PF cell ratio ' num2str(n_count_single/n_count_total)])
text(10,0.1, [ 'multi PF cell ratio ' num2str(n_count_multi/n_count_total)])

subplot(1,2,2)
hold on;
cdfplot(n_single_start_lap)
cdfplot(n_multi_start_lap)
title('Pf onset lap')
legend({'single','multi'})
[h p]=kstest2(n_single_start_lap,n_multi_start_lap);
text(10,0.2, [ 'kstest2 p=' num2str(p)])
%%
% instant PF single and multi PF ratio
A=sum(n_single_start_lap==1)/length(n_single_start_lap)
B=sum(n_multi_start_lap==1)/length(n_multi_start_lap)
%%
figure;
subplot(1,2,1);hold on;
bar(f_distribution)
grid off
title('f PF onset');

subplot(1,2,2)
n_distribution=histc(PF_start_lap,range);
hold on;bar(n_distribution)
ylim([1 180])
grid off


figure;
subplot(1,2,1)
hold on;
bar(f_distribution./sum(f_distribution))


grid off
subplot(1,2,2);
hold on;bar(n_distribution./sum(n_distribution))
 ylim([0 0.4])


norm_f=f_distribution/sum(f_distribution);
norm_n=n_distribution/sum(n_distribution);


%save('CA3_onset_dist','f_distribution','n_distribution','f_PF_start','n_PF_start');
% 
