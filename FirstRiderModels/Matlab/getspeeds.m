speeds = 0:1:10;
j = 1;
indices = [];
for i=1:length(v)
    if v(i) == speeds(j)
        indices = [indices;i];
        j = j+1;
    end
end
        