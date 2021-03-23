function mean_var=windowed_var(y,window)
total_var=[];
    for n=1:(length(y)-window+1)
        cur_var=var(y(:,n:n+window-1));
        total_var=[total_var cur_var];
    end
mean_var=var(total_var);


end