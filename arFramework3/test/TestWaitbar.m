
arWaitbar(0);
n = 1000;
delay = 0.01;

finaldelay = 0;
for j=1:n
    finaldelay = finaldelay + delay;
    delay = delay + 0.001;
end
disp(secToHMS(finaldelay))

delay = 0.01;
for j=1:n
    arWaitbar(j,n,'lala');
    pause(delay);
    delay = delay + 0.001;
end


arWaitbar(-1);
