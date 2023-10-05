function [d,v] = integrateGM(GM,g)

t = GM(:,1); dt = t(2)-t(1);
a = GM(:,2)*g;
v = zeros(length(t),1);
d = zeros(length(t),1);

% First integral
for i =1:length(t)-1
    v(i+1) = v(i) + (a(i+1)+a(i))*dt/2;
end

for i =1:length(t)-1
    d(i+1) = d(i) + (v(i+1)+v(i))*dt/2;
end

end