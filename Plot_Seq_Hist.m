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

bps_half_1=transpose(-N_half:-1);
bps_half_2=transpose(1:N_half);

for n=1:N_half
    bps(n,1)=bps_half_1(n,1);
end

for n=N_half+1:2*N_half
    bps(n,1)=bps_half_2(n-N_half,1);
end



x_line=-70:5:70;

for j=1:16;
    figure(j);
    h=plot(bps, P_seq(:,j),'r.-','LineWidth', 1.5, 'MarkerSize', 12);
    
    y_max=max(P_seq(:,j))+0.05;
    plot_lines(x_line,y_max)
    
    str_step=bps_string(j);
    title(str_step) 
    
    set(gca,'XTick',-70:10:70)
    ylim([0 y_max]);
    xlabel('Position (bp)')
    ylabel('Frequency')
    
    run make_string_fig;
    
    saveas(h,strfig_fig, 'fig')
    saveas(h,strfig_bmp, 'bmp')
end