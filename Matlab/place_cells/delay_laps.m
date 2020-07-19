clear all; %close all;
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
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)  
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
%             while start_lap==0
% %                 PF=PF_center(PF_count)-3:PF_center(PF_count)+3;
% %                 PF(PF<=0)=[];
%                 block=binmean_PF_temp(:,lap:lap+window-1);
%                 block(isnan(block))=0;
%                 
%                 if any(block(:,1)~=0)
%                 if sum(sum(block,1)~=0)<threshold
%                     lap=lap+1;
%                 else
%                     start_lap=lap;
%        
%                 end
%                 else
%                 lap=lap+1;
%                 end
%                 
%             end
            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
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




%% 
%clear all


%[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
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
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
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
%             while start_lap==0
% %                 PF=PF_center(PF_count)-3:PF_center(PF_count)+3;
% %                 PF(PF<=0)=[];
%                 block=binmean_PF_temp(:,lap:lap+window-1);
%                 block(isnan(block))=0;
%                 
%                 if any(block(:,1)~=0)
%                 if sum(sum(block,1)~=0)<threshold
%                     lap=lap+1;
%                 else
%                     start_lap=lap;
%        
%                 end
%                 else
%                 lap=lap+1;
%                 end
%                 
%             end
            PF_start_lap(PF_count)=start_lap;
            PF_count=PF_count+1;
        end
    
    end
                   
end
% n_PF_id=find(~isnan(number_of_PFs));
% instant_PF_id=n_PF_id(PF_start_lap==1);
% instantROI=cell_ROI(:,:,instant_PF_id);
% figure;
% imshow(reshape(mean(instantROI,3),[size(cell_ROI,1),size(cell_ROI,2)])*1.5);
% 
% delay_PF_id=n_PF_id(PF_start_lap>1);
% delayROI=cell_ROI(:,:,delay_PF_id);
% figure;
% imshow(reshape(mean(delayROI,3),[size(cell_ROI,1),size(cell_ROI,2)]).*2);

end
 display(['n day2 percent' ]);
 display(sum(n_total_firstlap_act>0)/length(n_total_firstlap_act));
 n_PF_start=PF_start_lap;
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

% figure;hold on;
% bar(f_distribution)
% b2=bar(n_distribution,'FaceColor',[0,0.7,0.7])
% set(get(b2,'Children'),'FaceAlpha',0.5)

figure;
subplot(1,2,1)
hold on;
bar(f_distribution./sum(f_distribution))

% n_distribution=histc(PF_start_lap,range);
% hold on;bar(n_distribution)
grid off
subplot(1,2,2);
hold on;bar(n_distribution./sum(n_distribution))
 ylim([0 0.4])
% grid off
% title('n PF onset');
%legend({'f','n'})
%ylim([0 350]);

%figure
norm_f=f_distribution/sum(f_distribution);
norm_n=n_distribution/sum(n_distribution);

% cdfplot(f_PF_start(f_PF_start>0));
% hold on; cdfplot(n_PF_start(n_PF_start>0));
% xlim([0 30])
% p=ranksum(f_PF_start(f_PF_start>0),n_PF_start(n_PF_start>0));
% title(['CA3 wilcoxon test,p= ' num2str(p)]);
% xlabel('PC onset lap')
% legend({'f','n'})
% grid off

save('CA3_onset_dist','f_distribution','n_distribution','f_PF_start','n_PF_start');
% 
% %%
% instant_PF_id=n_PF_id(PF_start_lap==1);
% instantROI=cell_ROI(:,:,instant_PF_id);
% figure;
% imshow(reshape(mean(instantROI,3),[size(cell_ROI,1),size(cell_ROI,2)])*1.5);
% %%
% instant_PF_id=n_PF_id(PF_start_lap>1);
% instantROI=cell_ROI(:,:,instant_PF_id);
% figure;
% imshow(reshape(mean(instantROI,3),[size(cell_ROI,1),size(cell_ROI,2)]).*2);