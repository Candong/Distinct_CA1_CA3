function [x y x_M y_M start_lap_v standard_lap_COM]=linreg_for_bylap_groupedCOM(f_sig_PFs,pf_id,norm_COMshift,startlap,stepsize,standardlap,stoplap)
% calculate windowed COM of PF acroos all animals 
% x:lap number y: pf_shift
y=[];
x=[];
jump_y=[];
jump_x=[];
x_M=[];
y_M=[];
start_lap_v=[];
standard_lap_COM=[];
 for i =1:size(f_sig_PFs,2)
     for p=1:size(f_sig_PFs,1)
       cur_sig_PFs=f_sig_PFs{p,i};
       cur_COMshift=norm_COMshift{p,i};
       use_id=find(cur_COMshift>-10& cur_COMshift<10); 
       %use_id=1:length(startlap{p,i});
       %lapnum=size(cur_sig_PFs{1,use_id(1)},2);  
       cur_pf_id=pf_id{p,i};
       for j=use_id
            cur_startlap=startlap{p,i}(j);
            if cur_startlap<20
            lapnum=size(cur_sig_PFs{1,cur_pf_id(j)},2);  
            lapnum=min(lapnum,stoplap+stepsize);
            cur_binmean=cur_sig_PFs{1,cur_pf_id(j)};
            if sum(mean(cur_binmean(:,1:5),2)>0)<51 && ~isclipped(cur_binmean)
            cur_COMv=NaN(1,stoplap+1);
            cur_jump_y=[];
            cur_jump_x=[];
            for n=1:(lapnum-stepsize+1)
                %if n>= cur_startlap
                step_binmean=cur_binmean(:,n:n+stepsize-1);
                [meanCOM SP]=meanCOMandSP(step_binmean',1,50,50);
                cur_COMv(n)=meanCOM;
%             if mod(n,stepsize)==1
%                cur_jump_y=[cur_jump_y meanCOM];
%                cur_jump_x=[cur_jump_x n];
%               
%                 
%             end
                %end
            end  
            cur_y=cur_COMv-cur_COMv(standardlap);
            %cur_y=diff(cur_COMv)
            cur_x=1:stoplap+1;
            y_M=[y_M; cur_y];
            x_M=[x_M; cur_x];
            start_lap_v=[start_lap_v cur_startlap];
            standard_lap_COM=[standard_lap_COM cur_COMv(standardlap)];
            remove_id=isnan(cur_y);
            cur_y(remove_id)=[];                    
            cur_x(remove_id)=[];
            y=[y cur_y];
            x=[x cur_x];
            
%             cur_jump_y=cur_jump_y-cur_COMv(standardlap+1);
%             jump_remove_id=isnan(cur_jump_y);
%             cur_jump_y(jump_remove_id)=[];
%             cur_jump_x(jump_remove_id)=[];
%             jump_y=[jump_y cur_jump_y];
%             jump_x=[jump_x cur_jump_x];
            end
            end    
        end
        
       end
         
     end
    
 end
