N=147;

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

bps=1:N-1;

for j=1:16;
    figure(j);
    h=plot(bps, P_seq(:,j),'r.-','LineWidth', 1.5, 'MarkerSize', 12);
    
    y_max=max(P_seq(:,j))+0.05;
    
    str_step=bps_string(j);
    title(str_step) 
    
    set(gca,'XTick',10:10:140)
    ylim([0 y_max]);
    xlabel('Position (bp)')
    ylabel('Frequency')
    
    run make_string_fig;
    
    %saveas(h,strfig_fig, 'fig')
    %saveas(h,strfig_bmp, 'bmp')
end