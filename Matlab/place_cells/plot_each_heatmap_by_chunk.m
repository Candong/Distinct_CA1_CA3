clear all;
close all;

[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
load([temp behavior_filepaths]);

[cellsort_filepaths, temp]=uigetfile('*.mat', 'Chose cellsort files to load:','MultiSelect','on');
load([temp cellsort_filepaths]);

%% change the frame numbers according to your recorded data
%analyze specific range:
startF=[1 5001 10001 14001 18001];
endF=[5000 10000 14000 18000 22000];
Fc3_DF=data.Fc3;
ybinned=behavior.ybinned;
ybinned=ybinned';
numbins=40;
cur_shock=behavior.shock(startF(end):endF(end));
shock_id=cur_shock>3;


%day1F_E=ca
for neuron_id=51:60
    figure('Position',[100,100,1200,400]);
    for i =1:length(startF)
        cur_E=calculate_E(ybinned(startF(i):endF(i)),0.01,0.04);
        binMean=chucked_heatmap(ybinned(startF(i):endF(i)),Fc3_DF(startF(i):endF(i),neuron_id),cur_E,numbins);
        
        subplot(1,length(startF),i)
        imagesc(binMean')
        if i ==length(startF)
            shock_lap=unique(cur_E(shock_id));
            shock_lap(shock_lap==0)=[];
            hold on;
            for j=shock_lap
            line([0 numbins],[j j],'Color','red')
            end
            
            
        end

    end 
end
%%


% function E=calculate_E(ybinned)
% ybinned_GoodBehav=ybinned;
% 
% %%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%
% 
% E=bwlabel(double_thresh(ybinned,0.036,0.035)); %labels each traversal
% % figure;plot(E)
% % hold on;plot(ybinned)
% trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
% trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units
% %ZThreshDivisor=4;%used to further break up the data in cases of long silent stretches
% %binsize=(trackend-trackstart)/numbins;%in Quake unit
% ybinmax=0.612;
% %% correct E
% wrong_lap=0;
%         for i=1:max(E)
% 
%         if i==1
%         onpoint=find(E==i,1);
%         offpoint=find(E==i,1,'last');  
%         if onpoint==1 & ybinned(onpoint)>0.12
%             E(onpoint:offpoint)=0;
%         end
% %         if ybinned(onpoint-1)>0.12
% %             E(E~=0)=E(E~=0)-1;
% %         end
%         
%         end
% 
%     end
% for i=1:max(E)
%     onpoint=find(E==i,1);
%     offpoint=find(E==i,1,'last');    
% %         end
%     count=1;
%     while max(ybinned(onpoint:offpoint))<ybinmax
%        count=count+1;
%        if i~=max(E) 
%          offpoint=find(E==i+1,1,'last');
%          E(onpoint:offpoint)=i;
%          E(offpoint+1:end)=E(offpoint+1:end)-1;
%        else
%            E(onpoint:offpoint)=0;
%        end
% 
%     end
% 
% end
% 
% 
% 
% 
% end

function binMean=chucked_heatmap(ybinned,Fc3_DF,E,numbins)
ybinned_GoodBehav=ybinned;
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units
%ZThreshDivisor=4;%used to further break up the data in cases of long silent stretches
binsize=(trackend-trackstart)/numbins;%in Quake unit

cut_mat=NaN(2000,1000);
cut_mat_ybinned=NaN(2000,1000);
wrong_lap=0;
for i=1:max(E)

    onpoint=find(E==i,1);
    offpoint=find(E==i,1,'last');    

    if max(ybinned(onpoint:offpoint))>0.6
    %cut out the transient
    cut_mat(1:(offpoint-onpoint)+1,i)=Fc3_DF(onpoint:offpoint); %cut out activity from each lap
    %cut out the ybinned associated with the transient
    cut_mat_ybinned(1:(offpoint-onpoint)+1,i)=ybinned(onpoint:offpoint,1);
    %figure; hold on;plot(cut_mat_ybinned);
    wrong_lap=0;
    else
    wrong_lap=1;    
    end

end
for i=1:max(E)
    [cut_mat_ybinned_sorted(:,i), idx] = sort(cut_mat_ybinned(:,i),1);
    cut_mat_sorted(:,i)=cut_mat(idx,i);
end

%bin activity
cut_mat_ybinned_sorted(cut_mat_ybinned_sorted==0)=nan;
topEdge = trackend; % define limits
botEdge = trackstart; % define limits
binMean=[];
for r=1:max(E)
    x = cut_mat_ybinned_sorted(:,r); %split into x and y
    y = cut_mat_sorted(:,r);
    binEdges = linspace(botEdge, topEdge, numbins+1);
    [h,whichBin] = histc(x, binEdges);
    %whichBin
    for i = 1:numbins
        flagBinMembers = (whichBin == i);
        binMembers     = y(flagBinMembers);
        binMean(i,r)     = mean(binMembers);
    end
end




end