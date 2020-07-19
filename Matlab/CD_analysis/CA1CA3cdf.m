
ca1_fpf_start=CA1_onset.f_PF_start;
ca1_npf_start=CA1_onset.n_PF_start;
ca3_fpf_start=CA3_onset.f_PF_start;
ca3_npf_start=CA3_onset.n_PF_start;

figure;
cdfplot(ca1_fpf_start(ca1_fpf_start>0));
hold on; cdfplot(ca3_fpf_start(ca3_fpf_start>0));
xlim([0 30])
p=ranksum(ca1_fpf_start(ca1_fpf_start>0),ca3_fpf_start(ca3_fpf_start>0));
title(['CA1 CA3 f, wt est,p= ' num2str(p)]);
xlabel('PC onset lap')
legend({'CA1','CA3'})
grid off

figure;
cdfplot(ca1_npf_start(ca1_npf_start>0));
hold on; cdfplot(ca3_npf_start(ca3_npf_start>0));
xlim([0 30])
p=ranksum(ca1_npf_start(ca1_npf_start>0),ca3_npf_start(ca3_npf_start>0));
title(['CA1 CA3 n, wt est,p= ' num2str(p)]);
xlabel('PC onset lap')
legend({'CA1','CA3'})
grid off

