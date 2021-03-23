function [f_delta_com f_lap]= shapecell2vector_COM_lap(f_all_deltaCOM)
f_delta_com=[];
f_lap=[];
for f=1:size(f_all_deltaCOM,2)
    cur_f_dcom=f_all_deltaCOM{f};
    cur_f_lap=ones(size(cur_f_dcom,1),size(cur_f_dcom,2)).*(1:size(cur_f_dcom,2));
    cur_f_dcom=reshape(cur_f_dcom,1,[]);
    cur_f_lap=reshape(cur_f_lap,1,[]);
    f_remove_id=isnan(cur_f_dcom);
    cur_f_dcom(f_remove_id)=[];
    cur_f_lap(f_remove_id)=[];
    f_delta_com=[f_delta_com cur_f_dcom];
    f_lap=[f_lap cur_f_lap];
    %cur_n_dcom=n_all_deltaCOM{f};
    
    
    
end
end