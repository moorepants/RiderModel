eigorg(:,1)=eigval(:,1); % set first column in organized matrix to the first column of the original matrix
unused=[1:10]; % create vector of unused indices
for i = 1:length(v)-1
    for j=1:10
        first = eigorg(j,i); % look at the first entry
        for k=1:length(unused)
            second = eigval(unused(k),i+1); % look at the entries in the next column
            diff(k) = abs(real(second) - real(first)); % calculate the difference in the previous column and the next column
        end
        [mindiff ind]=diff; % find smallest difference
        x=ind*ones(length(unused))-unused(m);
        for m=1:length(unused)
            if x == 0
                unused(m,1)=[]; % create 
            end
        end