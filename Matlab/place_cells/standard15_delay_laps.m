clear all; %close all;
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;
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
range=-1:100;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            %binmean_temp=sig_PFs_with_noise{j,i};
            binmean_PF_temp=sig_PFs{j,i};            
            mean_PF=mean(binmean_PF_temp,2);
            %mean_PF_withnoise=mean(binmean_temp,2);
            start_bin=PF_start_bins(j,i);
            end_bin=PF_end_bins(j,i);
            max_lap=size(binmean_PF_temp,2);
            
            
            start_lap=find_delaylap(binmean_PF_temp,start_bin,end_bin,window,threshold,max_lap);
            
            
           
            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
        end
    
    end
                   
end
end
f_distribution=histc(PF_start_lap,range);
figure;bar(f_distribution)
title('familiar');
%ylim([0 350]);

%% 
%clear all


[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;
for f=1:size(behavior,2)
    if contains(behavior{f},'_n_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);PF_center=[];
PF_center=[];
PF_mean=[];
PF_count=1;
% window=6;
% threshold=3;
PF_start_lap=[];
% range=-1:100;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            %binmean_temp=sig_PFs_with_noise{j,i};
            binmean_PF_temp=sig_PFs{j,i};            
            mean_PF=mean(binmean_PF_temp,2);
            %mean_PF_withnoise=mean(binmean_temp,2);
            start_bin=PF_start_bins(j,i);
            end_bin=PF_end_bins(j,i);
            max_lap=size(binmean_PF_temp,2);
            
            
            start_lap=find_delaylap(binmean_PF_temp,start_bin,end_bin,window,threshold,max_lap);
            
            
           
            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
        end
    
    end
                   
end
end

n_distribution=histc(PF_start_lap,range);
figure;bar(n_distribution)
title('novel');
%ylim([0 300]);
%%
figure
norm_f=f_distribution/sum(f_distribution);
norm_n=n_distribution/sum(n_distribution);
cdfplot(norm_f);
hold on; cdfplot(norm_n);
