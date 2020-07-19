clear all
%close all


%load file
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
load([temp behavior_filepaths]);

[cellsort_filepaths, temp]=uigetfile('*.mat', 'Chose cellsort files to load:','MultiSelect','on');
load([temp cellsort_filepaths]);

%if ~contains(cellsort_filepaths,['P' behavior_filepaths(end-8:end-5) '_' behavior_filepaths(end-4:end-4)])
%error('The behavior file and Cellsort files are not consistent')
%end

%cell_ROI=data.icaSegments;
%% change the frame numbers according to your recorded data
%analyze specific range:
last=0;
tracklength=300; %in cm
%num_of_lap=1:25;
% % f10min
 %startF=1;
%startF=13001;
  %endF=8000;
% % 
% %n1st10min
startF=behavior.startframe(end);
endF=length(behavior.ybinned);
%endF=14000;
% % 
% %n2nd10min
% startF=behavior.startframe+7000;
%endF=length(behavior.ybinned);
% 
%flast5min
% startF=3500;
% endF=7000;

% %nlast5min
% startF=behavior.startframe+3500;
% endF=behavior.startframe+7000;


Fc3_DF=data.Fc3(startF:endF,:);
Fc2=data.Fc(startF:endF,:);
ybinned=behavior.ybinned(startF:endF);
empty_id=[];
% empty_id=(9162:9464)+13000-startF;
% Fc3_DF(empty_id,:)=[];
% Fc2(empty_id,:)=[];
% ybinned(:,empty_id)=[];
% 3-1 day1 13910:14110 day2 13560:14290; -7000 day3 9162:9464+13000
%4-3 day1 16950,17110


%%
ybinmax=0.612;
label_high=0.015;
label_low=0.0145;



ybinned_GoodBehav=ybinned;
Fc3_DF_GoodBehav=Fc3_DF;



%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%

ybinned=ybinned';
E=bwlabel(double_thresh(ybinned,label_high,label_low)); %labels each traversal
%w=NaN(200,2);
sig_PFs=cell(8,size(Fc2,2));
sig_PFs_with_noise=cell(8,size(Fc2,2));
numneurons=size(Fc3_DF,2);
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units
ZThreshDivisor=4;%used to further break up the data in cases of long silent stretches
total_pos_randfields_rew=zeros(8,size(Fc2,2));
PF_PVALS=zeros(8,size(Fc2,2));
PF_width=NaN(8,size(Fc2,2));
PF_rate=NaN(8,size(Fc2,2));
PF_start_bins=NaN(8,size(Fc2,2));
PF_end_bins=NaN(8,size(Fc2,2));
number_of_PFs=NaN(1,size(Fc2,2));


numbins=50;
numrand=1000;%num iterations for shuffle
binsize=(trackend-trackstart)/numbins;%in Quake unit

width_M=zeros(1000,1000);

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
        if onpoint==1
            onpoint=2;
        end
        if ybinned(onpoint-1)>0.12
            E(E~=0)=E(E~=0)-1;
        end
        
        end

    end
%num_of_lap=(max(E)-25):max(E); last=1;   
% E(E<min(num_of_lap)|E>max(num_of_lap))=0;
% E=E-min(num_of_lap)+1;

 %E(E>num_of_lap)=0;
figure;plot(E)
%hold on;plot(behavior.ybinned)
hold on;plot(ybinned)


mean_trans=[];
binMean_M=[];

for ii = 1:numneurons
    
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
            if sum(flagBinMembers==1)>0
            binMembers     = y(flagBinMembers);
            binMean(i,r)     = mean(binMembers);
%             else
%                 binMean(i,r)=0;
            end
        end
    end
    
    figure;imagesc(binMean')
    
    mean_trans(ii,:)=mean(binMean,2);
    binMean_M(:,:,ii)=binMean';
    %remove possible errors from transients that occur during
    %teleportation:
%     for i=1:max(E)
%         errors=[];
%         %% might need to change
%         errors=bwlabel(double_thresh(binMean(:,i),0.0001,0.00001));%label a single run
%         for r=1:max(errors);
%             num=[];
%             num=find(errors==r);
%             if size(num,2)<=1 %transients that are this number or less are removed
%                 binMean(num,i)=0; %if tranient has 3 values or less, make zero
%             end
%         end
%     end
%     binMean(isnan(binMean))=0;
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    minfieldwidth=20; %in cm
    maxfieldwidth=150; %in cm
    minratio=3;%4;%original =3
    thresh1=0.15; %only coninous parts of the mean place field above this value are included
    minDF=0.1; %orignal = 0.1 %peak of the mean PF must be larger than this value
    minrate=15;
    %minrate=0.4;%0.4;%original =.3 %percentage of times a tranisent must occur in the PF
    baselength=0.25; %number of lowest values used to calculate a mean baseline value
    PVAL_thresh=0.05; %only PF that have pvals lower than this generated by shuffle test will be included
    
    parameters.minfieldwidth=minfieldwidth;
    parameters.maxfieldwidth=maxfieldwidth;
    parameters.minratio=minratio;
    parameters.thresh1=thresh1;
    parameters.minDF=minDF;
    parameters.minrate=minrate;
    parameters.baselength=baselength;
    parameters.PVAL_thresh=PVAL_thresh;
    
    parameters.startF=startF;
    parameters.endF=endF;
    
    
    %%search for place fields 
    
    %exclude non runing periods
    combinedYpos_rew=[ybinned Fc3_DF(:,ii)];
    sortedcombinedYpos_rew=sortrows(combinedYpos_rew);
    exclude_start=find(sortedcombinedYpos_rew(:,1)<trackstart);
    exclude_end=find(sortedcombinedYpos_rew(:,1)>trackend);
    exclude=[exclude_start;exclude_end];
    sortedcombinedYpos_rew(exclude,:)=[];
   
    
    %bin data:
    pos_Fpos_rew=[];
    x=sortedcombinedYpos_rew(:,1);
    y=sortedcombinedYpos_rew(:,2);
    [h,whichBin] = histc(x, binEdges);
    for j=1:numbins
        flagBinMembers = (whichBin == j);
        binMembers     = y(flagBinMembers);
        pos_Fpos_rew(j)     = mean(binMembers);
    end

    
    %tempypos_rew=pos_Fpos_rew(1,:); %contains the binned y position
    pos_rate_rew=[];
    fieldindpos_rew=zeros(numneurons,8);
    
    fieldcount=1;
    
    temp1=smooth(pos_Fpos_rew,3); %temp1 is the smoothed mean PF of the current neuron
    temp1sort=sort(temp1);
    base1=mean(temp1sort(1:round(baselength*numbins)));
    thresh2=(thresh1*(max(temp1)-mean(base1)))+base1;
    temp2=temp1>thresh2; %indices of values above threshold in the mean PF
    temp3=temp1<thresh2; %indices of out of PF firing
    L=bwlabeln(temp2,4);%labels each potential PF above the threshold
    num_potential_PFs=max(L);
    num_PFs=0;
    shuffle_num=0;
    
    %%%%%%%%%%%loops through each potential PF%%%%%%%%%%%%%%%%%%%%%%%%
    
    for tt=1:num_potential_PFs
        PF_start=[];
        PF_end=[];
        
        
        %calculate a new baseline and threshold for each potential PF, to account for differences in PF amplitudes
        
        thresh2=(thresh1*(max(temp1(find(L==tt)))-mean(base1)))+base1; %continuous region of PF has to be above this threshold (25% above the difference between the max bin and the mean of lowest 25 values - baselength)
        temp2=temp1>thresh2; %indices of values above threshold in the mean PF
        temp3=temp1<thresh2; %indices of out of PF firing
        
        
        
        if ((max(temp1(find(L==tt)))))>minDF
            %width=(((trackend-trackstart)/100)*length(find(L==tt)))*(150/(trackend-trackstart));%for 300cm track, converts to cm
            width=length(find(L==tt))*binsize*(tracklength/(trackend-trackstart));
            %width_M(tt,ii)=width;
            if width>minfieldwidth
            if width<maxfieldwidth    
                %rate=(size(binMean,2) - (sum(all(binMean(find(L==tt),:)==0,1)))) / size(binMean,2);
                rate=(size(binMean,2) - (sum(all(binMean(find(L==tt),:)==0,1))));
                width_M(tt,ii)=rate;
                if rate>minrate
                    if (mean(temp1(find(L==tt))))/mean(temp1(find(temp3==1)))>minratio %in PF firing has to be minratio above out of PF firing (excluding other PFs)
                        num_PFs=num_PFs+1; %counts number of potential PFs
                        PF_width(num_PFs,ii)=width; %in cm
                        PF_rate(num_PFs,ii)=rate;
                        PF_start=trackstart+((trackend-trackstart)/numbins)*min(find(L==tt));
                        PF_end=trackstart+((trackend-trackstart)/numbins)*max(find(L==tt));
                        %plot average image of cells
%                         image1=zeros(N,M,3);
%                         image1(:,:,1)=f0/max(max(f0));
%                         image1(:,:,2)=f0/max(max(f0));
%                         image1(:,:,3)=f0/max(max(f0));
%                         mask=squeeze(masks(ii,:,:));
%                         mask(mask >0) = 1;
%                         f=figure;image(image1);
%                         movegui(f,'northeast');
%                         hold on;contour(mask,'c','LineWidth',1);
%                         [row,col]=find(mask,1,'first');
%                         ROI = sprintf('%d', ii);
%                         text(col+7,row,ROI,'Color','r','FontWeight','Bold','FontSize',7);
                        %plot mean PF
%                         figure;
%                         hold on;
%                         plot(temp1,'LineWidth', 1.5);title(num2str(ii));
                        PF_start_bins(num_PFs,ii)=(PF_start - trackstart)/binsize;
                        PF_end_bins(num_PFs,ii)=(PF_end - trackstart)/binsize;
%                         plot(PF_start_bins(num_PFs,ii):PF_end_bins(num_PFs,ii),temp1(PF_start_bins(num_PFs,ii):PF_end_bins(num_PFs,ii)),'r','LineWidth', 1.5);
%                         %plot transients
                        %f=figure;imagesc(binMean');
%                         movegui(f,'northwest');
                        %pause (0.1);
                    end
                end
            end 
            end
        end
        
        if and(num_PFs==1, shuffle_num==0)
            %checks to see if there are any potential PFs. if yes, run shuffle.
            %if no, move to next PF. if more than 1, use PVAL from first PF
            shuffle_num=1;
            
            %%%%%%%%%%%%%original shuffle test + jims shuffle%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %take a PF one at a time, do a shuffle and test significance
            
            
            
            %%%%%%%%%%%%%%%%%%%jims code%%%%%%%%%%%%%%%%%%%%%%
            
            LLL = bwlabel(Fc3_DF(:,ii)>0);
            if max(LLL)>1
                for jj=1:max(LLL)-1
                    indsITI(jj) = find(LLL==jj+1,1,'first')-find(LLL==jj,1,'last');
                end
                zerosthresh = mean(indsITI);
            end
            
            if max(LLL)<=1
                zerosthresh = length(1:find(LLL==1,1,'first'));
            end
            
            %shuffle test:
            for hh=1:numrand
                %hh
                
                combinedYpos_rew=[];
                
                %randomly shuffle Fc3
                new_Fc3_DF=[];
                
                m = bwlabel(Fc3_DF(:,ii)>0);%this generates a vector M that relates the indicies from Fc3_DF_GoodBehav(:,ii) to the number of each sequential fluorescent trasient eg M=[000111000222200333
                counter = 1;
                MM = [1];
                
                for i = 1:length(m)-1%this generates a vector MM that use M and  Fc3_DF_GoodBehav(:,ii) eg MM=[111222333444455666];
                    if m(i)==m(i+1)
                        MM(i+1)=counter;
                    end
                    if m(i)~=m(i+1)
                        counter=counter+1;
                        MM(i+1)=counter;
                    end
                end
                
                P = MM;
                %%this section removes inter transient intervals with long ITIs
                for i = 1:max(MM)
                    ITI=length(find(MM==i));
                    if sum(Fc3_DF(find(MM==i),ii))==0 && ITI>zerosthresh/ZThreshDivisor
                        x=randi([1 ITI],1,1);
                        MM(find(MM==i,1,'first')+x:length(MM)) = MM(find(MM==i,1,'first')+x:length(MM))+1;%P(find(MM==i,1,'first')+x:length(MM)) = P(find(MM==i,1,'first')+x:length(MM))+1;
                        i=i-1;
                    end
                end
                
                x = randperm(max(MM));%randperm(max(P));x is a vector containing a randomly assorted set of the number of intra and inter transient intervals. eg. given MM, x could be x = [3 4 6 2 5 1]
                new_Fc3_DF=[];
                
                for i = 1:length(x)%this generates the shuffled data set new_Fc3_DF that containins a random order of each inter and intra transient intervals from Fc3_DF_GoodBehav(:,ii)
                    new_Fc3_DF = [new_Fc3_DF Fc3_DF(find(MM==x(i)),ii)'];%[new_Fc3_DF Fc3_DF_GoodBehav(find(P==x(i)),ii)'];
                end
                
                
                new_Fc3_DF=new_Fc3_DF'; %this is the shuffled trace (a new one is created with each iteration)
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%Dan's code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %bin position and average Ftrace
                combinedYpos_rew=[ybinned new_Fc3_DF];
                
                sortedcombinedYpos_rew=sortrows(combinedYpos_rew);
                exclude_start=find(sortedcombinedYpos_rew(:,1)<trackstart);
                exclude_end=find(sortedcombinedYpos_rew(:,1)>trackend);
                exclude=[exclude_start;exclude_end];
                sortedcombinedYpos_rew(exclude,:)=[];
                
                
                
   
                
                
                
                
                     %bin data shuffled data:
                     pos_Fpos_rew_shuffle=[];
                     x=sortedcombinedYpos_rew(:,1);
                     y=sortedcombinedYpos_rew(:,2);
                     [h,whichBin] = histc(x, binEdges);
                     for j=1:numbins
                         flagBinMembers = (whichBin == j);
                         binMembers     = y(flagBinMembers);
                         pos_Fpos_rew_shuffle(j)     = mean(binMembers);
                         numsigs_pos_rew(j)=sum(binMembers>0)/length(binMembers);
                     end
                
                
                %%search for place fields on the new shuffled trace
                fieldindpos_rew=zeros(numneurons,8);
                tempypos_rew_shuffle=pos_Fpos_rew_shuffle; %binned location
                fieldcount=1;
                temp1_shuffle=smooth(pos_Fpos_rew_shuffle,3); %binned dF of shuffled data
                temp1sort=sort(temp1_shuffle);
                base1=mean(temp1sort(1:round(baselength*numbins)));
                thresh2=(thresh1*(max(temp1_shuffle)-mean(base1)))+base1;
                temp2_shuffle=temp1_shuffle>thresh2;
                L_shuffle=bwlabeln(temp2_shuffle,4);
                
                
                for p=1:max(L_shuffle)
                    if mean(numsigs_pos_rew(find(L_shuffle==p)))>minrate
                        width=length(find(L==tt))*binsize*(tracklength/(trackend-trackstart));
                        %width=(((trackend-trackstart)/100)*length(find(L_shuffle==p)))*(tracklength/(trackend-trackstart));%for 300cm track, converts to cm
                        if ((mean(temp1(find(L_shuffle==p))))/mean(temp1(find(L_shuffle~=p))))>minratio
                            fieldindpos_rew(ii,fieldcount)=min(find(L_shuffle==p));
                            fieldindpos_rew(ii,fieldcount+1)=max(find(L_shuffle==p));
                            fieldcount=fieldcount+2;
                            if ((max(temp1(find(L_shuffle==p)))))>minDF
                                if width>minfieldwidth
                                if width<maxfieldwidth
                                    total_pos_randfields_rew(1,ii)=total_pos_randfields_rew(1,ii)+1;
%                                     figure;
%                                     hold on;
%                                     plot(temp1_shuffle,'LineWidth', 1.5); %displays shuffled data
%                                     ginput(1);
%                                     close all;
                                end
                                end
                            end
                        end
                    end
                end
                
                
            end %end of shuffle iteration cycle
            
            
            
        end %jumps to here if there is no potential PF (num_PFs will either be zero, 1 or more than 1.
        
        
        if num_PFs>0
            PF_PVALS(num_PFs,ii)=total_pos_randfields_rew(1,ii)/numrand; %each PF from the same neuron gets the same PVAL
            PF_PVALS(num_PFs+1:end,ii)=NaN; %make all othe elements equal NaN
            number_of_PFs(1,ii)=num_PFs;
        else
            PF_PVALS(:,ii)=NaN;%if there are no PFs make PVAL NaN
        end
        
        %if PVAL of the potential PF is avove Pval thresh, make PF features zero:
        if num_PFs>0
            if PF_PVALS(num_PFs,ii)>PVAL_thresh
                PF_width(num_PFs,ii)=0;
                PF_rate(num_PFs,ii)=0;
                PF_start_bins(num_PFs,ii)=0;
                PF_end_bins(num_PFs,ii)=0;
                number_of_PFs(1,ii)=0;
            end
        end
        
        %%%%%%%%%%Take significant PFs and zero out all transients that aren't part of the PF%%%%%%%%%
        if (num_PFs>0)
            if isempty(sig_PFs{num_PFs,ii})
                if PF_PVALS(num_PFs,ii)<PVAL_thresh
                    %cut out all transients that contribute to the mean significant PF
                    binmean_temp=binMean;
                    n=NaN(size(binmean_temp,1),size(binmean_temp,2));
                    for i=1:max(E)% loop through each lap
                        n(:,i)=bwlabeln(binmean_temp(:,i),4);
                        if any(n(PF_start_bins(num_PFs,ii):PF_end_bins(num_PFs,ii),i));%are any elements nonzero within PF limits?
                            %loop through each transient and make zeros
                            %if out of PF or keep if in PF
                            for h=1:max(n(:,i))
                                if any (find((n(PF_start_bins(num_PFs,ii):PF_end_bins(num_PFs,ii),i)==h)))    ;
                                    continue
                                else
                                    binmean_temp(find(n(:,i)==h),i)=0;
                                end
                            end
                        else
                            %make entire run zero
                            binmean_temp(:,i)=0;
                        end
                    end
                    %display significant PF
                    binmean_temp(isnan(binmean_temp))=0;
                    f=figure;imagesc(binmean_temp');
                    title(num2str(ii));
                    movegui(f,'southwest');
%                     f=figure;plot(smooth(mean(binmean_temp'),3),'r','LineWidth', 1.5);
%                     movegui(f,'south');
%                     pause (5); %pause to admire the significant PF
%                     ginput(1);
%                     close all;
                    %save significant PF
                    sig_PFs{num_PFs,ii}=binmean_temp;
                    sig_PFs_with_noise{num_PFs,ii}=binMean;
                end
            end
        end
        
        
    end %end of loop for cycling  through multiple PFs in the same neuron
    %close all
end %loops through each neuronn

% if startF>300
% if num_of_lap(1)==1
%  save([behavior_filepaths(1:end-4) '_n_PF_strict_1st25lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% elseif num_of_lap(1)==11
%     save([behavior_filepaths(1:end-4) '_n_PF_strict_late15lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% elseif num_of_lap(1)==26
%     save([behavior_filepaths(1:end-4) '_n_PF_strict_2nd25lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% end
% else
%     
% if last==1
%  save([behavior_filepaths(1:end-4) '_f_PF_strict_last25lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% elseif num_of_lap(1)==1
%  save([behavior_filepaths(1:end-4) '_f_PF_strict_1st25lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% elseif num_of_lap(1)==11
%     save([behavior_filepaths(1:end-4) '_f_PF_strict_late15lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% elseif num_of_lap(1)==26
%     save([behavior_filepaths(1:end-4) '_f_PF_strict_2nd25lap_10lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','empty_id','parameters');%,'-append');
% 
% end
% end
if startF>300
    save([behavior_filepaths(1:end-4) '_n_PF_strict_alllap_' num2str(minrate) 'lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','empty_id','parameters');%,'cell_ROI','-append');
else
    save([behavior_filepaths(1:end-4) '_f_PF_strict_alllap_' num2str(minrate) 'lap'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','empty_id','parameters');%,'cell_ROI','-append');
end

