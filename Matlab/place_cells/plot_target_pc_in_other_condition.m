clear all; close all;
[pf_filepaths, temp]=uigetfile('*.mat', 'Chose target files to load:','MultiSelect','on');
load([temp pf_filepaths])
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
load([temp behavior_filepaths])
[activity_filepaths, temp]=uigetfile('*.mat', 'Chose fc3 files to load:','MultiSelect','on');
load([temp activity_filepaths])

subplot_x=10;
subplot_y=ceil(sum(number_of_PFs(~isnan(number_of_PFs)))/subplot_x);
PC_count=0;
PC_id=find(~isnan(number_of_PFs));

startF=input(['What is the start frame number? ']);
endF=input(['What is the end frame number? ']);
ybinned=behavior.ybinned(startF:endF);
cur_E=calculate_E(ybinned,0.01,0.04);
cur_shock=behavior.airpuff(startF:endF);
%cur_shock=behavior.shock(startF:endF);
shock_id=cur_shock>3;

for i=PC_id%(1:100)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs{j,i};
            PC_count=PC_count+1;
        end
subplot(subplot_y,subplot_x,PC_count);
imagesc(binmean_temp');
title(['neuron id=' num2str(i)])

if any(shock_id)
shock_lap=unique(cur_E(shock_id));
shock_lap(shock_lap==0)=[];
hold on;
    for s=shock_lap
    line([0 size(binmean_temp,1)],[s s],'Color','red')
    end
end
    end
                
end
%%
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
load([temp behavior_filepaths])
[activity_filepaths, temp]=uigetfile('*.mat', 'Chose fc3 files to load:','MultiSelect','on');
load([temp activity_filepaths])

subplot_x=10;%10;
subplot_y=10;%ceil(sum(number_of_PFs(~isnan(number_of_PFs)))/subplot_x);

startF=input(['What is the start frame number? ']);
endF=input(['What is the end frame number? ']);

ybinned=behavior.ybinned(startF:endF);
E=calculate_E(ybinned,0.01,0.04);
Fc3_DF=data.Fc3(144001:end,:);
%Fc3_DF=data.Fc3(startF:endF,:);
cur_shock=behavior.airpuff(startF:endF);
%cur_shock=behavior.shock(startF:endF);
shock_id=cur_shock>3;
ybinned=behavior.ybinned(startF:endF);
ybinned=ybinned';
numbins=40;

trackstart=min(ybinned)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned)-0.005; %track end location in quake units
%ZThreshDivisor=4;%used to further break up the data in cases of long silent stretche
binsize=(trackend-trackstart)/numbins;%in Quake unit

%binsize=(trackend-trackstart)/numbins;%in Quake unit

width_M=zeros(1000,1000);

figure;
%%
count=1;
for ii =PC_id(1:100)%numneurons

    
    %seperate transients into laps
    cut_mat=NaN(2000,1000);
    cut_mat_ybinned=NaN(2000,1000);
    wrong_lap=0;
    for i=1:max(E)

        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');    

        if max(ybinned(onpoint:offpoint))>0.6
        %cut out the transient
        cut_mat(1:(offpoint-onpoint)+1,i)=Fc3_DF(onpoint:offpoint,ii); %cut out activity from each lap
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
        for i = 1:numbins
            flagBinMembers = (whichBin == i);
            binMembers     = y(flagBinMembers);
            binMean(i,r)     = mean(binMembers);
        end
    end

subplot(subplot_y,subplot_x,count);
imagesc(binMean');
title(['neuron id=' num2str(ii)])
count=count+1;
if any(shock_id)
shock_lap=unique(E(shock_id));
shock_lap(shock_lap==0)=[];
hold on;
    for s=shock_lap
    line([0 size(binmean_temp,1)],[s s],'Color','red')
    end
end
    

    %close all
end %loops through each neuron
