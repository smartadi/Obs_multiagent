function U = u_signal(L,seed)

th = 0.4;
avg = 5;
U=[];
for i = 1:L
    rng(seed+i)
    u = double(rand(1)>th);
    U = [U,u];
    if i <= avg
        U(:,end) = 0;
    elseif U(:,end-1) == 0 && mean(U(:,end-4:end)) > 1/avg
        U(:,end) = 0;
    end

end
