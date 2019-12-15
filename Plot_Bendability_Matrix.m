N=147;
N_half=73;

bond_min=1;
bond_max=14;

f1=0;
f2=0;
f3=0;

tail1=0;
tail2=0;

beta=3;
beta_seq=3;

N_Run=200;

sum_hist=zeros(N-1,16);

for Run=0:N_Run-1
    run make_string;
    hist=load(str_h_Seq, '-ascii');
    sum_hist=sum_hist+hist;
end

num=0;

for j=1:16
    num=num+sum_hist(1,j);
end

P_seq = sum_hist/num;

step_names={'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA','TC', 'TG', 'TT'};

imagesc(P_seq);
colorbar
shading faceted 
set(gca,'ytick',0:10:140)
set(gca,'xtick',1:16)
set(gca,'yMinorTick','on')
set(gca,'XTickLabel',step_names)
xlabel('Base pair Steps')
ylabel('Position (bp)')

%saveas(h,strfig_fig, 'fig')
%saveas(h,strfig_bmp, 'bmp')