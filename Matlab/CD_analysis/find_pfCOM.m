function COM_lap=find_pfCOM(binmean,stepsize)
%calculate by lap COM by stepsize
            maxlap=size(binmean,2);
            COM_lap=[];
            for i=1:maxlap-stepsize+1
                cur_lap_act =binmean(:,i:i+stepsize-1);
                [meanCOM SP]=meanCOMandSP(cur_lap_act',1,50,50);

                COM_lap=[COM_lap meanCOM];
                
      
            end


end