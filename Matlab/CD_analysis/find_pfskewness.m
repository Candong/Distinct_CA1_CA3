function skewness_lap=find_pfskewness(binmean,stepsize)
            maxlap=size(binmean,2);
            skewness_lap=[];
            for i=1:maxlap-stepsize+1
                block=binmean(:,i:i+stepsize-1);
                block_mean=mean(block,2);
                if any(block_mean>0)
                    cur_skewness=[];
                    for i=1:stepsize
                        cur_lap=block(:,i);
                        lap_transient = find(block(:,i)>0);
                        if ~isempty(lap_transient)
                        cur_skewness=[cur_skewness cal_skewness(cur_lap(lap_transient(1):lap_transient(end)))];
                        end

                    end 
                    cur_skewness=mean(cur_skewness);
                    
                else
                    cur_skewness=NaN;
                end
%                 if stepsize==1
%                     for i=1:maxlap
%                         block=binmean(:,i);
%                         if any(block>0)
%                         lap_transient = find(block>0);
%                         cur_skewness=cal_skewness(block(lap_transient(1):lap_transient(end)));
%                         else
%                             cur_skewness =NaN;
%                         end
%                         skewness_lap=[skewness_lap cur_skewness];
%                         
%                     end
%                     
%                 end
%                 cur_lap_act =mean(binmean(:,i:i+stepsize-1),2); % old
%                 %%%%%%%version
%                 block_transient=find(cur_lap_act>0);
%                 if ~isempty(block_transient)
%                     %cur_skewness=skewness(cur_lap_act(block_transient(1):block_transient(end)));
%                     cur_skewness=cal_skewness(cur_lap_act);
%                 else
%                     cur_skewness=NaN;
%                 end
% 
                skewness_lap=[skewness_lap cur_skewness];
                
      
            end


end

