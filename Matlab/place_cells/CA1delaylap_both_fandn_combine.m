clear all; close all;

% delay lap of familiar PF found by last 15 laps
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
PF_info=struct;
field_name={'plane1','plane2','plane3'};

PF_mean=[];
PF_count=1;
PF_start_lap=[];

planenum=2;

for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
place_cell_id=[];
for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            place_cell_id=[place_cell_id i];
        end
    end   
end


PF_info(f).filepaths=behavior_filepaths{f};
PF_info(f).place_cell_id=place_cell_id; 
PF_info(f).PF_start_bin=PF_start_bins;
PF_info(f).PF_end_bin=PF_end_bins;
PF_info(f).number_of_PFs=number_of_PFs;

end
%% parameters
window=6;
threshold=3;
%% delap laps in f
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_late15_pfid{1}=PF_info(i).place_cell_id;
          f_late15_id{1}=i;
          f_late15_pfnum{1}=PF_info(i).number_of_PFs;
%           elseif contains(PF_info(i).filepaths,'plain2')
%           f_late15_pfid{2}=PF_info(i).place_cell_id;
%           f_late15_id{2}=i;
%           f_late15_pfnum{2}=PF_info(i).number_of_PFs;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_late15_pfid{2}=PF_info(i).place_cell_id;
          f_late15_id{2}=i;
          f_late15_pfnum{2}=PF_info(i).number_of_PFs;
          end
          
      elseif contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_1st25_pfid{1}=PF_info(i).place_cell_id;
          f_1st25_id{1}=i;
          f_1st25_pfnum{1}=PF_info(i).number_of_PFs;
%           elseif contains(PF_info(i).filepaths,'plain2')
%           f_1st25_pfid{2}=PF_info(i).place_cell_id;
%           f_1st25_id{2}=i;
%           f_1st25_pfnum{2}=PF_info(i).number_of_PFs;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_1st25_pfid{2}=PF_info(i).place_cell_id;
          f_1st25_id{2}=i;      
          f_1st25_pfnum{2}=PF_info(i).number_of_PFs;
          end
          
      end
      
  end
    
end

%% find pcell from first 25 or late 15
for i=1:planenum
f_potential_samepf{i}=intersect(f_1st25_pfid{i},f_late15_pfid{i});
f_1st25multifield{i}=find(f_1st25_pfnum{i}>1);
f_late15multifield{i}=find(f_late15_pfnum{i}>1);

end
% make sure all potential same pf is the same
for p=1:planenum
for i=f_potential_samepf{p} 
    if ~isempty(i)
% load([temp behavior_filepaths{f_1st25_id{p}}]);

if (~ismember(i,f_1st25multifield{p}&~ismember(i,f_late15multifield{p})))
   
       f_1st25_pf(1)=PF_info(f_1st25_id{p}).PF_start_bin(1,i);
       f_1st25_pf(2)=PF_info(f_1st25_id{p}).PF_end_bin(1,i); 

% load([temp behavior_filepaths{f_late15_id{p}}]);    
     
       f_late15_pf(1)=PF_info(f_late15_id{p}).PF_start_bin(1,i);
       f_late15_pf(2)=PF_info(f_late15_id{p}).PF_end_bin(1,i);

pf_diff=sum(abs(f_1st25_pf-f_late15_pf));

if pf_diff>6
    remove_pf=find(f_potential_samepf{p}==i);
    f_potential_samepf{p}(remove_pf(1))=[];
end

else 
   pf_diff=[]; 
   for j=1:2
       f_1st25_pf(1)=PF_info(f_1st25_id{p}).PF_start_bin(j,i);
       f_1st25_pf(2)=PF_info(f_1st25_id{p}).PF_end_bin(j,i); 

% load([temp behavior_filepaths{f_late15_id{p}}]);    
     
       f_late15_pf(1)=PF_info(f_late15_id{p}).PF_start_bin(j,i);
       f_late15_pf(2)=PF_info(f_late15_id{p}).PF_end_bin(j,i);  
       
   pf_diff(j)=sum(abs(f_1st25_pf-f_late15_pf));
   
   if pf_diff(j)>6
    remove_pf=find(f_potential_samepf{p}==i);
    f_potential_samepf{p}(remove_pf(j))=[];
   end
end  
    end
    end
end
PF_info(f_1st25_id{p}).same_pfid=f_potential_samepf{p};

end
%% find delay lap for f
% window=6;
% threshold=2;
range=-1:50;
figure;
for plane=1:planenum
 behaviorfile=dir(['*ds_plain' num2str(plane) '.mat']);
 load(behaviorfile.name);
 cellsortfile=dir(['*Plain_' num2str(plane) '_*MotCor*' '.mat']);
 load(cellsortfile.name);
 neuron_id=combine2array(f_1st25_pfid{plane},f_potential_samepf{plane}, f_late15_pfid{plane});
f_binmean=[];
 %find start and end_bin for each place cell
 start_bin=[];end_bin=[];
 first25_id=reshape(setxor(f_1st25_pfid{plane},f_potential_samepf{plane}),1,[]);
 late15_id=f_late15_pfid{plane};
 if ~isempty(first25_id)
      load([temp behavior_filepaths{f_1st25_id{plane}}]);    
 for k=unique(first25_id)

    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(f_1st25_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(f_1st25_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
 end
  end
%  
   load([temp behavior_filepaths{f_late15_id{plane}}]);
  for k=unique(late15_id)

    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(f_late15_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(f_late15_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
  end
 
 PF_info(f_late15_id{plane}).combinePC=neuron_id;
 PF_info(f_late15_id{plane}).combinePC_start=start_bin;
 PF_info(f_late15_id{plane}).combinePC_end=end_bin;
 PF_info(f_late15_id{plane}).PC_ROI=cell_ROI(:,:,neuron_id);
 
 %analyze specific range:
 num_of_lap=1:25;
startF=behavior.startframe(1);
startF=1;
endF=13000;
Fc3_DF=data.Fc3(startF:endF,:);
Fc2=data.Fc(startF:endF,:);
ybinned=behavior.ybinned(startF:endF);
% 
empty_id=[];
% empty_id=(13610:14290)-startF;
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
        if ybinned(onpoint)>0.11
            E(E~=0)=E(E~=0)-1;
        end
        
        end

    end
E(E<min(num_of_lap)|E>max(num_of_lap))=0;
E=E-min(num_of_lap)+1;
mean_trans=[];
%%
count=1;
for ii =neuron_id
    
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
%    figure;
%    imagesc(binMean');
%    title(num2str(ii));

    delay_lapid{plane}(count)=find_delaylap(binMean,start_bin(count),end_bin(count),window,threshold,max(E));
  denoise_binmean=binMean';
  denoise_binmean(:,round([1:start_bin(count) end_bin(count): numbins]))=0;
  
  f_binmean(:,:,count)=denoise_binmean;
count=count+1;
    %close all
end %loops through each neuron
 
PF_info(f_1st25_id{plane}).delay_lap=delay_lapid{plane};
PF_info(f_1st25_id{plane}).multipf_neuron=f_1st25multifield{plane};
PF_info(f_late15_id{plane}).multipf_neuron=f_late15multifield{plane};
PF_info(f_1st25_id{plane}).combinebinMean=f_binmean;

f_distribution=histc(delay_lapid{plane},range);
PF_info(f_1st25_id{plane}).delay_lap_distribute=f_distribution;

if ~isempty(f_distribution)
subplot(1,3,plane)
bar(f_distribution)
title(['familiar plane' num2str(plane)]);
%ylim([0 25]);
end
 
end



%% nlast 15 and 1st 25
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_n_')
      
      if contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_1st25_pfid{1}=PF_info(i).place_cell_id;
          n_1st25_id{1}=i;
          n_1st25_pfnum{1}=PF_info(i).number_of_PFs;
%           elseif contains(PF_info(i).filepaths,'plain2')
%           n_1st25_pfid{2}=PF_info(i).place_cell_id;
%           n_1st25_id{2}=i;          
%           n_1st25_pfnum{2}=PF_info(i).number_of_PFs;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_1st25_pfid{2}=PF_info(i).place_cell_id;
          n_1st25_id{2}=i; 
          n_1st25_pfnum{2}=PF_info(i).number_of_PFs;
          end
          
      elseif contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_late15_pfid{1}=PF_info(i).place_cell_id;
          n_late15_id{1}=i;
          n_late15_pfnum{1}=PF_info(i).number_of_PFs;
%           elseif contains(PF_info(i).filepaths,'plain2')
%           n_late15_pfid{2}=PF_info(i).place_cell_id;
%           n_late15_id{2}=i;
%           n_late15_pfnum{2}=PF_info(i).number_of_PFs;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_late15_pfid{2}=PF_info(i).place_cell_id;
          n_late15_id{2}=i;
          n_late15_pfnum{2}=PF_info(i).number_of_PFs;
          end
          
      end
      
  end
    
end

for i=1:planenum
n_potential_samepf{i}=intersect(n_1st25_pfid{i},n_late15_pfid{i});
n_1st25multifield{i}=find(n_1st25_pfnum{i}>1);
n_late15multifield{i}=find(n_late15_pfnum{i}>1);
end
% make sure all potential same pf is the same
for p=1:planenum
for i=n_potential_samepf{p} 
    if ~isempty(i)
% load([temp behavior_filepaths{f_1st25_id{p}}]);

if (~ismember(i,n_1st25multifield{p}&~ismember(i,n_late15multifield{p})))
   
       n_1st25_pf(1)=PF_info(n_1st25_id{p}).PF_start_bin(1,i);
       n_1st25_pf(2)=PF_info(n_1st25_id{p}).PF_end_bin(1,i); 
 
     
       n_late15_pf(1)=PF_info(n_late15_id{p}).PF_start_bin(1,i);
       n_late15_pf(2)=PF_info(n_late15_id{p}).PF_end_bin(1,i);

pf_diff=sum(abs(n_1st25_pf-n_late15_pf));

if pf_diff>6
    remove_pf=find(n_potential_samepf{p}==i);
    n_potential_samepf{p}(remove_pf(1))=[];
end

else 
   pf_diff=[]; 
   for j=1:2
       n_1st25_pf(1)=PF_info(n_1st25_id{p}).PF_start_bin(j,i);
       n_1st25_pf(2)=PF_info(n_1st25_id{p}).PF_end_bin(j,i); 

     
       n_late15_pf(1)=PF_info(n_late15_id{p}).PF_start_bin(j,i);
       n_late15_pf(2)=PF_info(n_late15_id{p}).PF_end_bin(j,i);  
       
   pf_diff(j)=sum(abs(n_1st25_pf-n_late15_pf));
   
   if pf_diff(j)>6
    remove_pf=find(n_potential_samepf{p}==i);
    n_potential_samepf{p}(remove_pf(j))=[];
   end
end  
    end
    end
end
PF_info(n_1st25_id{p}).same_pfid=n_potential_samepf{p};

end

%% find delay laps
% window=6;
% threshold=2;
delay_lapid=cell(1,3);
range=-1:50;
ii=[];
figure;
for plane=1:planenum
 behaviorfile=dir(['*ds_plain' num2str(plane) '.mat']);
 load(behaviorfile.name);
 cellsortfile=dir(['*Plain_' num2str(plane) '_*MotCor*' '.mat']);
 load(cellsortfile.name);
 neuron_id=combine2array(n_1st25_pfid{plane},n_potential_samepf{plane}, n_late15_pfid{plane});

 %find start and end_bin for each place cell
 start_bin=[];end_bin=[];
 first25_id=reshape(setxor(n_1st25_pfid{plane},n_potential_samepf{plane}),1,[]);
 late15_id=n_late15_pfid{plane};
 n_binmean=[];
 if ~isempty(first25_id)
  load([temp behavior_filepaths{n_1st25_id{plane}}]);
 for k=unique(first25_id)

    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(n_1st25_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(n_1st25_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
 end
 end
 
 if ~isempty(late15_id)
load([temp behavior_filepaths{n_late15_id{plane}}]);
  for k=unique(late15_id)

    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(n_late15_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(n_late15_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
  end
 end
 PF_info(n_late15_id{plane}).combinePC=neuron_id;
 PF_info(n_late15_id{plane}).combinePC_start=start_bin;
 PF_info(n_late15_id{plane}).combinePC_end=end_bin;
 PF_info(n_late15_id{plane}).PC_ROI=cell_ROI(:,:,neuron_id);
 
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
% empty_id=(13590:14280)-startF;
% Fc3_DF(empty_id,:)=[];
% Fc2(empty_id,:)=[];
% ybinned(:,empty_id)=[];

%3-1 day1 13910:14110 day2 13590:14280; 4-3 day1 16950,17110

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
for ii =neuron_id
    
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
%    figure;
%    imagesc(binMean');
%    title(num2str(ii));

    delay_lapid{plane}(count)=find_delaylap(binMean,start_bin(count),end_bin(count),window,threshold,max(E));
    denoise_binmean=binMean';
    denoise_binmean(:,round([1:start_bin(count) end_bin(count): numbins]))=0;
  
    n_binmean(:,:,count)=denoise_binmean;
count=count+1;
    %close all
end %loops through each neuron
 
PF_info(n_1st25_id{plane}).delay_lap=delay_lapid{plane};
PF_info(n_1st25_id{plane}).multipf_neuron=n_1st25multifield{plane};
PF_info(n_late15_id{plane}).multipf_neuron=n_late15multifield{plane};
PF_info(n_1st25_id{plane}).combinebinMean=n_binmean;
 
n_distribution=histc(delay_lapid{plane},range);
PF_info(n_1st25_id{plane}).delay_lap_distribute=n_distribution;

if ~isempty(n_distribution)
subplot(1,3,plane)
bar(n_distribution)
title(['novel plane' num2str(plane)]);
%ylim([0 450]);
end
 
end
PF_info(1).empty_id=empty_id;
save([PF_info(1).filepaths(1:17) '_info_threshold=' num2str(threshold)],'PF_info','-v7.3')

%% function for change negetive num to 0

function a=rectl(num)
if num<=0
    a=0;
else
    a=num;
end


end


function newarray=combine2array(A,B,D)

C=setxor(A,B);
C=A(ismember(A,C));

if size(C,1)>1
    C=C';
end
newarray=[C D];
if isempty(newarray)
    newarray=[];
end
end
