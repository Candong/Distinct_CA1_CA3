
function width_lap=find_pfwidth(binmean,window)
            maxlap=size(binmean,2);
            width_lap=[];lap=1;
                while lap+window-1<=maxlap
                block=binmean(:,lap:lap+window-1);
                block_mean=mean(block,2);
                if any(block_mean>0)
                    cur_width=[];
                    for i=1:window
                        lap_transient = find(block(:,i)>0);
                        if ~isempty(lap_transient)
                        cur_width=[cur_width lap_transient(end)-lap_transient(1)+1];
                        end

                    end 
                    cur_width=mean(cur_width);
%                     block_transient=find(block_mean>0); % old way to
%                     cur_width= block_transient(end)-block_transient(1)+1;
                else
                    cur_width=0;
                end
                width_lap=[width_lap cur_width];
                
                lap=lap+1;
                end
%                 if length(width_lap)>=25 & width_lap (25 ~=0) & width_lap(1)~=0
%                 if width_lap (25)-width_lap(1)>10
%                     figure;
%                     plot(mean(binmean(:,1:5),2));
%                     hold on;
%                     plot(mean(binmean(:,25:29),2));
%                 end
 
%                end

end