function [x y  x_M y_M ]=linreg_for_bylap_groupedCOM_onsetlap(f_sig_PFs,pf_id,norm_COMshift,startlap,stepsize,standardlap,stoplap)
% calculate windowed COM of PF acroos all animals 
% x:lap number y: pf_shift
y=[];
x=[];
jump_y=[];
jump_x=[];
x_M=[];
y_M=[];
% x_all_lap=[];
% y_all_lap=[];
 for i =1:size(f_sig_PFs,2)
     for p=1:size(f_sig_PFs,1)
       cur_sig_PFs=f_sig_PFs{p,i};
       cur_COMshift=norm_COMshift{p,i};
       use_id=find(cur_COMshift>-10& cur_COMshift<10);  
       %lapnum=size(cur_sig_PFs{1,use_id(1)},2);  
       cur_pf_id=pf_id{p,i};
       for j=use_id
            cur_startlap=startlap{p,i}(j);
            if cur_startlap>0
            
            %lapnum=min(lapnum,stoplap+stepsize);
            cur_binmean=cur_sig_PFs{1,cur_pf_id(j)};
            cur_binmean=cur_binmean(:,cur_startlap:end);
            lapnum=size(cur_binmean,2);  
            lapnum=min(lapnum,stoplap+stepsize);
            if  ~isclipped(cur_binmean) & lapnum+1>stepsize & cur_startlap>1%sum(mean(cur_binmean(:,1:5),2)>0)<51 &&
            
            cur_COMv=NaN(1,stoplap+1);

              
            for n=1:(lapnum-stepsize+1)
                step_binmean=cur_binmean(:,n:n+stepsize-1);
                [meanCOM SP]=meanCOMandSP(step_binmean',1,50,50);
                cur_COMv(n)=meanCOM;
            end
            
            
            
            cur_y=cur_COMv-cur_COMv(standardlap);
            cur_x=1:stoplap+1;
            y_M=[y_M; cur_y];
            x_M=[x_M; cur_x];
            remove_id=isnan(cur_y);
            cur_y(remove_id)=[];                    
            cur_x(remove_id)=[];
            y=[y cur_y];
            x=[x cur_x];

            end
            end    
        end
     end 
       end
         
     end

