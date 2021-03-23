function E=calculate_E(ybinned,range_1,range_2)
ybinned_GoodBehav=ybinned;

%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%

E=bwlabel(double_thresh(ybinned,range_1,range_2)); %labels each traversal
% figure;plot(E)
% hold on;plot(ybinned)
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units
%ZThreshDivisor=4;%used to further break up the data in cases of long silent stretches
%binsize=(trackend-trackstart)/numbins;%in Quake unit
ybinmax=0.61;
%% correct E
wrong_lap=0;
        for i=1:max(E)

        if i==1
        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');  
        if onpoint==1 & ybinned(onpoint)>0.12
            E(onpoint:offpoint)=0;
        end
%         if ybinned(onpoint-1)>0.12
%             E(E~=0)=E(E~=0)-1;
%         end
        
        end

    end
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




end