function state=double_thresh(data,thresh1,thresh2)
%state=double_thresh(data,thresh1,thresh2)
%thresh1 upperthreshold
%thresh2 lowerthreshold

state=zeros(1,length(data));
if data(1)<thresh2
    state(1)=0;
else
    state(1)=1;
end

for i=2:length(data)
    state(i)=state(i-1);
    if state(i-1)==0
        if data(i)>thresh1
            state(i)=1;
        end
    else
        if data(i)<thresh2
            state(i)=0;
        end
    end
end
