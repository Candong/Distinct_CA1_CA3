clear all; %close all;
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
% put in all CA1 or CA3 place field files, no matter f or n, those files
% with the' ... 15_lap' end.
file_count=1;

if ~isa(behavior,'cell') 
  behavior={behavior};  
    
end 

%% Familiar
for f=1:size(behavior,2)
    if contains(behavior{f},'_f_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);



PF_mean=[];
f_total_SD=[];
f_clip_id=[];
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean=sig_PFs{j,i};
            SD=cal_meanPF_SD(mean(binmean,2));
            f_total_SD=[f_total_SD SD]; % calculate the SD
            
            if isclipped(binmean)
                f_clip_id=[f_clip_id  1]; % find out if it is a clipped PF or not'
            else
                f_clip_id=[f_clip_id 0];
            end


        end
    
    end
                   
end

end




%% Novel
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
n_total_SD=[];
n_clip_id=[];
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean=sig_PFs{j,i};
            SD=cal_meanPF_SD(mean(binmean,2));
            n_total_SD=[n_total_SD SD];
            
            if isclipped(binmean)
                n_clip_id=[n_clip_id  1];
            else
                n_clip_id=[n_clip_id 0];
            end


        end
    
    end
                   
end

end
%%
figure;
histogram(n_total_SD,0:1:120,'Normalization','Probability')
line([median(n_total_SD) median(n_total_SD)], [0,0.03],'Color','r')
title('CA1 n')
text(100,0.03,num2str(median(n_total_SD)))

n_remove_clip_SD=n_total_SD(~n_clip_id);
figure;
histogram(n_remove_clip_SD,0:1:120,'Normalization','Probability')
line([median(n_remove_clip_SD) median(n_remove_clip_SD)], [0,0.03],'Color','r')
text(100,0.03,num2str(median(n_remove_clip_SD)))
ylim([0 0.09])
title('CA1 n no-clipped')