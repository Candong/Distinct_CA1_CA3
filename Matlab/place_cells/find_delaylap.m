function start_lap=find_delaylap(binmean,start_bin,end_bin,window,threshold,max_lap)
            start_lap=0;lap=1;
            while start_lap==0 %& lap+window-1<max_lap
                if lap+window-1<max_lap
                block=binmean(start_bin:end_bin,lap:lap+window-1);
                block(isnan(block))=0;
                
                if any(block(:,1)~=0)
                if sum(sum(block,1)~=0)<threshold
                    lap=lap+1;
                else
                    start_lap=lap;
       
                end
                else
                lap=lap+1;
                end
                else
                    start_lap=-1;
                end
                
            end

end