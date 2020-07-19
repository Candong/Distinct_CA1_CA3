
% function[e_PF_sorted]=lay_PF_by_lap(lap_num )
clear all;
%prompt = 'What is the lap numbers you want to use to plot this map? ';
[filepaths, temp]=sort(uigetfile('*.mat', 'Chose n files to load:','MultiSelect','on'));
e_lap_num = 1:20;
m_lap_num =30:50;
l_lap_num =60:80;
e_PF_center=[];
e_PF_mean=[];
m_PF_mean=[];
l_PF_mean=[];
PF_count=1;
placecell_id={};
inconsistent_id=[];
MEAN=[];
if ~isa(filepaths,'cell') 
  filepaths={filepaths};  
    
end 

for f=1:size(filepaths,2)
load([filepaths{f}]);
placecell_id{f}=[];
for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs_with_noise{j,i};
            binmean_temp(isnan(binmean_temp))=0;
            binmean_PF_temp=sig_PFs{j,i};
            binmean_PF_temp(isnan(binmean_PF_temp))=0;
            
            e_mean_PF=mean(binmean_PF_temp(:,e_lap_num),2);
            e_mean_PF_withnoise=mean(binmean_temp(:,e_lap_num),2);
            
            m_mean_PF=mean(binmean_PF_temp(:,m_lap_num),2);
            m_mean_PF_withnoise=mean(binmean_temp(:,m_lap_num),2);
            l_mean_PF=mean(binmean_PF_temp(:,l_lap_num),2);
            l_mean_PF_withnoise=mean(binmean_temp(:,l_lap_num),2);
            
            if (sum(e_mean_PF)~=0 & sum(m_mean_PF)~=0 & sum(l_mean_PF))
            [e_PF_center(PF_count) SP]=meanCOMandSP(binmean_PF_temp(:,e_lap_num)',1,50,50);
            [m_PF_center(PF_count) SP]=meanCOMandSP(binmean_PF_temp(:,m_lap_num)',1,50,50);
            [l_PF_center(PF_count) SP]=meanCOMandSP(binmean_PF_temp(:,l_lap_num)',1,50,50);
            %[meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
            e_PF_mean(PF_count,:)=e_mean_PF_withnoise;
            m_PF_mean(PF_count,:)=m_mean_PF_withnoise;
            l_PF_mean(PF_count,:)=l_mean_PF_withnoise;
            
            PF_count=PF_count+1;

            else

            end
            
        end
    
    end
                   
end
%MEAN=[MEAN; mean_trans(placecell_id{f},:)];
end
%%
[e_PF_sorted,e_sort_id]=sort(e_PF_center);
e_PF_mean_sorted=e_PF_mean(e_sort_id,:);

remove_id=find(mean(e_PF_mean_sorted')>3);
e_PF_mean_sorted(remove_id,:)=[];

[m_PF_sorted,m_sort_id]=sort(m_PF_center);
m_PF_mean_sorted=m_PF_mean(e_sort_id,:);

remove_id=find(mean(m_PF_mean_sorted')>3);
m_PF_mean_sorted(remove_id,:)=[];

[l_PF_sorted,l_sort_id]=sort(l_PF_center);
l_PF_mean_sorted=l_PF_mean(e_sort_id,:);

remove_id=find(mean(l_PF_mean_sorted')>3);
l_PF_mean_sorted(remove_id,:)=[];

%%
f=figure;
subplot(1,3,1)
imagesc((e_PF_mean_sorted'./max(e_PF_mean_sorted'))');
title(['sorted PFs form lap ' num2str(e_lap_num(1)) 'to lap ' num2str(e_lap_num(end))]);
colormap jet

subplot(1,3,2)
imagesc((m_PF_mean_sorted'./max(m_PF_mean_sorted'))');
title(['sorted PFs form lap ' num2str(m_lap_num(1)) 'to lap ' num2str(m_lap_num(end))]);
colormap jet

subplot(1,3,3)
imagesc((l_PF_mean_sorted'./max(l_PF_mean_sorted'))');
title(['sorted PFs form lap ' num2str(l_lap_num(1)) 'to lap ' num2str(l_lap_num(end))]);
colormap jet
%%
figure;plot(flip(e_PF_sorted),1:length(e_PF_center))
hold on;
plot(flip(m_PF_sorted),1:length(e_PF_center))
plot(flip(l_PF_sorted),1:length(e_PF_center))
legend({'early','middle','late'})

%%
figure;plot((e_PF_sorted),(1:length(e_PF_center)),'.','MarkerSize',10)
hold on;
plot((m_PF_sorted),(1:length(e_PF_center)),'.','MarkerSize',10)
plot((l_PF_sorted),(1:length(e_PF_center)),'.','MarkerSize',10)
ylim([0 length(e_PF_center)])
legend({'early','middle','late'})
set(gca, 'YDir','reverse')
%%
% m_PF_center=[];
% m_PF_mean=[];
% PF_count=1;
% %placecell_id={};
% MEAN=[];
% m_lap_num = 30:50;
% for f=1:size(filepaths,2)
% load([filepaths{f}]);
% %placecell_id{f}=[];
% cur_pc_id=placecell_id{f};
% for i=cur_pc_id
%     for j=1:size(sig_PFs,1)
%         if ~isempty(sig_PFs{j,i})
%             binmean_temp=sig_PFs_with_noise{j,i};
%             binmean_temp(isnan(binmean_temp))=0;
%             binmean_PF_temp=sig_PFs{j,i};
%             binmean_PF_temp(isnan(binmean_PF_temp))=0;
%             
%             mean_PF=mean(binmean_PF_temp(:,m_lap_num),2);
%             mean_PF_withnoise=mean(binmean_temp(:,m_lap_num),2);
%             if sum(mean_PF)~=0
%             [m_PF_center(PF_count) SP]=meanCOMandSP(binmean_PF_temp(:,m_lap_num)',1,50,50);
%             %[meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
%             m_PF_mean(PF_count,:)=mean_PF_withnoise;
%             PF_count=PF_count+1;
%             %placecell_id{f}=[placecell_id{f} i];
%             end
%         end
%     
%     end
%                    
% end
% %MEAN=[MEAN; mean_trans(placecell_id{f},:)];
% end
% 
% [m_PF_sorted,sort_id]=sort(m_PF_center);
% m_PF_mean_sorted=m_PF_mean(sort_id,:);
% familiar_sort=sort_id;
% remove_id=find(mean(m_PF_mean_sorted')>3);
% m_PF_mean_sorted(remove_id,:)=[];
% f=figure;imagesc(m_PF_mean_sorted);
% title(['sorted PFs form lap ' num2str(m_lap_num(1)) 'to lap ' num2str(m_lap_num(end))]);
% colormap jet







%%
% l_PF_center=[];
% l_PF_mean=[];
% PF_count=1;
% %placecell_id={};
% MEAN=[];
% l_lap_num = 60:80;
% for f=1:size(filepaths,2)
% load([filepaths{f}]);
% %placecell_id{f}=[];
% cur_pc_id=placecell_id{f};
% for i=cur_pc_id
%     for j=1:size(sig_PFs,1)
%         if ~isempty(sig_PFs{j,i})
%             binmean_temp=sig_PFs_with_noise{j,i};
%             binmean_temp(isnan(binmean_temp))=0;
%             binmean_PF_temp=sig_PFs{j,i};
%             binmean_PF_temp(isnan(binmean_PF_temp))=0;
%             
%             mean_PF=mean(binmean_PF_temp(:,l_lap_num),2);
%             mean_PF_withnoise=mean(binmean_temp(:,l_lap_num),2);
%             if sum(mean_PF)~=0
%             [l_PF_center(PF_count) SP]=meanCOMandSP(binmean_PF_temp(:,l_lap_num)',1,50,50);
%             %[meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
%             l_PF_mean(PF_count,:)=mean_PF_withnoise;
%             PF_count=PF_count+1;
%             %placecell_id{f}=[placecell_id{f} i];
%             end
%         end
%     
%     end
%                    
% end
% %MEAN=[MEAN; mean_trans(placecell_id{f},:)];
% end
% 
% [l_PF_sorted,sort_id]=sort(l_PF_center);
% l_PF_mean_sorted=l_PF_mean(sort_id,:);
% familiar_sort=sort_id;
% remove_id=find(mean(l_PF_mean_sorted')>3);
% l_PF_mean_sorted(remove_id,:)=[];
% f=figure;imagesc(l_PF_mean_sorted);
% title(['sorted PFs form lap ' num2str(l_lap_num(1)) 'to lap ' num2str(l_lap_num(end))]);
% colormap jet


% %[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose n files to load:','MultiSelect','on');
% PF_center=[];
% PF_mean=[];
% PF_count=1;
% % if ~isa(behavior_filepaths,'cell') 
% %   behavior_filepaths={behavior_filepaths};  
% %     
% % end 
% for f=1:size(n_filepaths,2)
% load([temp n_filepaths{f}]);
% for i=1:size(sig_PFs,2)
%     for j=1:size(sig_PFs,1)
%         if ~isempty(sig_PFs{j,i})
%             binmean_temp=sig_PFs_with_noise{j,i};
%             binmean_PF_temp=sig_PFs{j,i};
%             
%             mean_PF=mean(binmean_PF_temp,2);
%             mean_PF_withnoise=mean(binmean_temp,2);
%             
%             PF_center(PF_count)=find(mean_PF==max(mean_PF));
%             PF_mean(PF_count,:)=mean_PF_withnoise;
%             PF_count=PF_count+1;
% 
%         end
%     
%     end
%                    
% end
% end
% [PF_sorted,sort_id]=sort(PF_center);
% PF_mean_sorted=PF_mean(sort_id,:);
% 
% f=figure;imagesc(PF_mean_sorted);
% title('novel day2')
% colormap jet
% % f=figure;imagesc(PF_mean_sorted./max(PF_mean_sorted'));
% % title('novel')
% 
% MEAN=[];
% 
% for f=1:size(n_filepaths,2)
% load([temp n_filepaths{f}]);
% 
% field_id=placecell_id{f};
% MEAN=[MEAN; mean_trans(placecell_id{f},:)];
% 
% end
% f2n_mean=MEAN(familiar_sort,:);
% f2n_mean(remove_id,:)=[];
% f=figure;imagesc(f2n_mean)%,[0 3.4]);
% title('day1 to day2')
% colormap jet

%end