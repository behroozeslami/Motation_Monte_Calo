function plot_lines(x,y_max)
hold on
for j=1:length(x)
    x_array=[x(j) x(j)];
    y_array=[0 y_max];
    
    str_marker='k--';
    if round(j/2)==j/2
        str_marker='k-';
    end
    
    plot(x_array,y_array,str_marker, 'LineWidth', 1.5)
end
hold off