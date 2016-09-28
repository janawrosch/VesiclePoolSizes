function [ x_werte ] = scale_x_axis( n_pics )

x_werte=zeros(1,n_pics);
count=0;
for i=1:28
    count=count+1;
    x_werte(1,i)=count;
    
end
for i=29:42
    count=count+0.25;
    x_werte(1,i)=count;
end
for i=43:82
    count=count+1;
    x_werte(1,i)=count;
end
for i=83:197
    count=count+0.25;
    x_werte(1,i)=count;
end
for i=198:220
    count=count+1;
    x_werte(1,i)=count;
end
for i=221:306
    count=count+0.25;
    x_werte(1,i)=count;
end
for i=307:n_pics
    count=count+1;
    x_werte(1,i)=count;
end

x_werte2=zeros(1,n_pics);
count=0;
for i=1:28
    count=count+1;
    x_werte2(1,i)=count;
    
end
for i=29:42
    count=count+0.25;
    x_werte2(1,i)=count;
end
for i=43:82
    count=count+1;
    x_werte2(1,i)=count;
end
for i=83:197
    count=count+0.25;
    x_werte2(1,i)=count;
end
for i=198:270
    count=count+0.3;
    x_werte2(1,i)=count;
end
for i=271:306
    count=count+0.25;
    x_werte2(1,i)=count;
end
for i=307:n_pics
    count=count+1;
    x_werte2(1,i)=count;
end


end

