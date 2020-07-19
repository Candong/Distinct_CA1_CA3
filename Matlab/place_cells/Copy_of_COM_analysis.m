% this script calculate the COM of each PF
clear all; close all;
%COM familiar

% delay lap of familiar PF found by last 15 laps
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose info files to load:','MultiSelect','on');
%PF_info=struct;
field_name={'plane1','plane2','plane3'};

PF_mean=[];
PF_count=1;

PF_start_lap=[];


load([temp behavior_filepaths]);

%% parameters
window=6;
threshold=2;

%% delap laps in f
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_late15_pfid{1}=PF_info(i).combinePC;
          f_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_late15_pfid{2}=PF_info(i).combinePC;
          f_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          f_late15_pfid{3}=PF_info(i).combinePC;
          f_late15_id{3}=i;      
          end
          
      end
      
  end
    
end

%% find delay COM for f

for plane=1:3
 behaviorfile=dir(['*ds_plain' num2str(plane) '.mat']);
 load(behaviorfile.name);
 cellsortfile=dir(['*Plain_' num2str(plane) '_*MotCor*' '.mat']);
 load(cellsortfile.name);
 %neuron_id=combine2array(n_first25_pfid{plane},n_potential_samepf{plane}, n_late15_pfid{plane});

 %find start and end_bin for each place cell
 start_bin=[];end_bin=[];

 if ~isempty(f_late15_pfid{plane})

 start_bin=PF_info(f_late15_id{plane}).PF_start_bin;
 end_bin=PF_info(f_late15_id{plane}).PF_end_bin;
 
 end
 
 
 %PF_info(f_late15_id{plane}).combinePC=neuron_id;
%  PF_info(f_late15_id{plane}).combinePC_start=start_bin;
%  PF_info(f_late15_id{plane}).combinePC_end=end_bin;
 
 %analyze specific range:
 num_of_lap=1:25;
% startF=behavior.startframe(1);
% endF=length(behavior.ybinned);
startF=1;
endF=7000;
Fc3_DF=data.Fc3(startF:endF,:);
Fc2=data.Fc(startF:endF,:);
ybinned=behavior.ybinned(startF:endF);
% 
% empty_id=13610:14290-startF;
% Fc3_DF(empty_id,:)=[];
% Fc2(empty_id,:)=[];
% ybinned(empty_id,:)=[];



%%

ybinmax=0.612;
label_high=0.012;
label_low=0.011;


ybinned_GoodBehav=ybinned;
Fc3_DF_GoodBehav=Fc3_DF;



%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%

ybinned=ybinned';
E=bwlabel(double_thresh(ybinned,label_high,label_low)); %labels each traversal

%numneurons=size(Fc3_DF,2);
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units

numbins=50;

binsize=(trackend-trackstart)/numbins;%in Quake unit


%% correct E
    wrong_lap=0;
    for i=1:max(E)

            
        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');    
%         end
        count=1;
        while max(ybinned(onpoint:offpoint))<ybinmax 
           count=count+1;
           if i~=max(E) 
             offpoint=find(E==i+1,1,'last');
             E(onpoint:offpoint)=i;
             E(offpoint+1:end)=E(offpoint+1:end)-1;
           else
               E(onpoint:offpoint)=0;
           end

        end
        
    end

    
    for i=1:max(E)

        if i==1
        onpoint=find(E==i,1);
        if ybinned(onpoint)>0.08
            E(E~=0)=E(E~=0)-1;
        end
        
        end

    end
E(E<min(num_of_lap)|E>max(num_of_lap))=0;
E=E-min(num_of_lap)+1;
mean_trans=[];
%%
count=1;
pf_COM=[];
pf_SP=[];
for ii = f_late15_pfid{plane}
    
    %%%%%%%%now cut out transients from each lap and bin them%%%%%%%%%%%%%%%%%
    
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
    

    figure;
   imagesc(binMean');
   title([num2str(plane) ' ' num2str(ii)])
%%%%%%%%%%%%%%%%%%%%%%
PF_withoutnoise=binMean';

for pfs=sum(~isnan(start_bin(:,ii)))
start_b=start_bin(pfs,ii);
end_b=end_bin(pfs,ii);
[meanCOM cur_SP]=meanCOMandSP(PF_withoutnoise,start_b,end_b,numbins);
pf_COM=[meanCOM pf_COM]; 
pf_SP=[cur_SP pf_SP];
%find(sum(PF_withoutnoise(:,start_b:end_b),2)==0)
PF_empty{pfs,ii}=find(sum(PF_withoutnoise(:,start_b:end_b),2)==0);
PF_denoise=PF_withoutnoise;
PF_denoise(:,round([1:start_b-1 end_b-1:numbins]))=0;
figure;imagesc(PF_denoise');
sig_PF{pfs,ii}=PF_denoise;
end
f_PF_empty{plane}=PF_empty;
f_sig_PF{plane}=sig_PF;
end
clear PF_empty
clear sig_PF
f_COM{plane}=pf_COM;
f_SP{plane}=pf_SP;


end
%PF_withoutnoise(:,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COM in n
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_n_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_late15_pfid{1}=PF_info(i).place_cell_id;
          n_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_late15_pfid{2}=PF_info(i).place_cell_id;
          n_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_late15_pfid{3}=PF_info(i).place_cell_id;
          n_late15_id{3}=i;      
          end
          
      end
      
  end
    
end

%% find delay lap for n
n_delay_lapid=cell(1,3);
% window=6;
% threshold=2;
% range=-1:50;
% figure;
for plane=1:3
 behaviorfile=dir(['*ds_plain' num2str(plane) '.mat']);
 load(behaviorfile.name);
 cellsortfile=dir(['*Plain_' num2str(plane) '_*MotCor*' '.mat']);
 load(cellsortfile.name);
 %neuron_id=combine2array(n_first25_pfid{plane},n_potential_samepf{plane}, n_late15_pfid{plane});

 %find start and end_bin for each place cell

 start_bin=[];end_bin=[];

 if ~isempty(n_late15_pfid{plane})

 start_bin=PF_info(n_late15_id{plane}).PF_start_bin;
 end_bin=PF_info(n_late15_id{plane}).PF_end_bin;
 
 end

 
 
 %PF_info(f_late15_id{plane}).combinePC=neuron_id;
%  PF_info(n_late15_id{plane}).combinePC_start=start_bin;
%  PF_info(n_late15_id{plane}).combinePC_end=end_bin;
 
 %analyze specific range:
 num_of_lap=1:25;
startF=behavior.startframe(end);
endF=length(behavior.ybinned);
% startF=1;
% endF=7000;
Fc3_DF=data.Fc3(startF:endF,:);
Fc2=data.Fc(startF:endF,:);
ybinned=behavior.ybinned(startF:endF);
% 
empty_id=[];
% empty_id=(16950:17110)-startF;
% Fc3_DF(empty_id,:)=[];
% Fc2(empty_id,:)=[];
% ybinned(:,empty_id)=[];


%%

ybinmax=0.612;
label_high=0.012;
label_low=0.011;


ybinned_GoodBehav=ybinned;
Fc3_DF_GoodBehav=Fc3_DF;



%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%

ybinned=ybinned';
E=bwlabel(double_thresh(ybinned,label_high,label_low)); %labels each traversal

%numneurons=size(Fc3_DF,2);
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units

numbins=50;

binsize=(trackend-trackstart)/numbins;%in Quake unit


%% correct E
    wrong_lap=0;
    for i=1:max(E)

            
        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');    
%         end
        count=1;
        while max(ybinned(onpoint:offpoint))<0.638 
           count=count+1;
           if i~=max(E) 
             offpoint=find(E==i+1,1,'last');
             E(onpoint:offpoint)=i;
             E(offpoint+1:end)=E(offpoint+1:end)-1;
           else
               E(onpoint:offpoint)=0;
           end

        end
        
    end

    
    for i=1:max(E)

        if i==1
        onpoint=find(E==i,1);
        if ybinned(onpoint)>0.02
            E(E~=0)=E(E~=0)-1;
        end
        
        end

    end
E(E<min(num_of_lap)|E>max(num_of_lap))=0;
E=E-min(num_of_lap)+1;
mean_trans=[];
%%
count=1;
pf_COM=[];
pf_SP=[];
for ii = n_late15_pfid{plane}
    
    %%%%%%%%now cut out transients from each lap and bin them%%%%%%%%%%%%%%%%%
    
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
    
%     mean_trans(ii,:)=mean(binMean,2);
%     figure;
%     imagesc(binMean');
%     title('n');
PF_withoutnoise=binMean';

for pfs=sum(~isnan(start_bin(:,ii)))
start_b=start_bin(pfs,ii);
end_b=end_bin(pfs,ii);
[meanCOM cur_SP]=meanCOMandSP(PF_withoutnoise,start_b,end_b,numbins);
pf_COM=[meanCOM pf_COM]; 
pf_SP=[cur_SP pf_SP];
%find(sum(PF_withoutnoise(:,start_b:end_b))==0)
PF_empty{pfs,ii}=find(sum(PF_withoutnoise(:,start_b:end_b),2)==0);
PF_denoise=PF_withoutnoise;
PF_denoise(:,round([1:start_b-1 end_b-1:numbins]))=0;
sig_PF{pfs,ii}=PF_denoise;
end
n_PF_empty{plane}=PF_empty;
n_sig_PF{plane}=sig_PF;

end
clear PF_empty
clear sig_PF
n_COM{plane}=pf_COM;
n_SP{plane}=pf_SP;

%     n_delay_lapid{plane}(count)=find_delaylap(binMean,start_bin(count),end_bin(count),window,threshold,max(num_of_lap));
%     
% count=count+1;
    %close all
end %loops through each neuron
 save([PF_info(1).filepaths(1:17) '_COM_and_SP'],'f_SP','f_COM','n_SP','f_COM','f_late15_pfid','n_late15_pfid','f_PF_empty','n_PF_empty','f_sig_PF','n_sig_PF');
%     PF_info(n_late15_id{plane}).delay_lap=n_delay_lapid{plane};
%  
% n_distribution=histc(n_delay_lapid{plane},range);
% PF_info(n_late15_id{plane}).delay_lap_distribute=n_distribution;
% 
% if ~isempty(n_distribution)
% subplot(1,3,plane)
% bar(n_distribution)
% title(['novel plane' num2str(plane)]);
% ylim([0 25]);
% end
%  
 
% end
 


