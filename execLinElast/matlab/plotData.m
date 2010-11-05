plot(data.x,data.tau/4,'.-','Color',C(1,:))
hold on
plot(data.x,data.V,'.-','Color',C(2,:))
hold off
colors
hold on
for i = 1:max(data.L)
    I = data.L == i;
    plot(data.x(I),i+0*data.tau(I)/4,'.','Color',C(i,:))
   %  plot(data.x(I),data.V(I),'.','Color',C(i,:))
end
hold off