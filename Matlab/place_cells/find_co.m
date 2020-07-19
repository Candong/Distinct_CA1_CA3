function c=find_co(id1,id2,num_pf1,num_pf2,sig_PFs)
% num_pf1=1;
% num_pf2=1;
% if id1==id2
%     num_pf2=2;
% end

a=mean(sig_PFs{num_pf1,id1},2);
b=mean(sig_PFs{num_pf2,id2},2);
c=corr(a,b);


end