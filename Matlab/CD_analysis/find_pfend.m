function end_lap=find_pfend(binmean,stepsize,find_start)
            maxlap=size(binmean,2);
            end_lap=[];
            for i=1:maxlap-stepsize+1
                cur_lap_act =mean(binmean(:,i:i+stepsize-1),2);
                block_transient=find(cur_lap_act>0);
                if ~isempty(block_transient)
                    if find_start==0
                    cur_end=block_transient(end);
                    else
                    cur_end=block_transient(1);
                    end
                else
                    cur_end=nan;
                end

                end_lap=[end_lap cur_end];
                
      
            end


end