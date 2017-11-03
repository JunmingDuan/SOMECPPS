hold on
load data;
%xlabel('log(xn/x1)')
%ylabel('log(yn/y1)')
x1=0.01:0.01:0.99;
%x2=-1:0.02:1;
plot(x1,data,'ob');
%plot(x1,data(2,1:201),'-g')
%plot(x2,data(3,1:101),'*k')
%plot(x2,data(4,1:101),'*m')
%plot(x1,data(5,1:201),'--b')
%plot(x1,data(6,1:201),'--g')
%plot(x2,data(7,1:101),'+k')
%plot(x2,data(8,1:101),'+m')
plot(x1,0.5*(1-exp(-x1/1))/(1-exp(-1/1))+0.5*x1,'-r');
%plot(x1,exp(-0.5*((x1-3).^2)/(sqrt(1+2*0.001*0.15)))/sqrt(1+2*0.001*0.15),'r')
%legend('h=0.005,mu=0.1,EX','h=0.005,mu=0.5,EX','h=0.01,mu=0.1,EX','h=0.01,mu=0.5,EX',...
    %'h=0.005,mu=0.1,UW','h=0.005,mu=0.5,UW','h=0.01,mu=0.1,UW','h=0.01,mu=0.5,UW','exact')
%legend('delta=e-3,t=0.1,EX','delta=e-3,t=0.5,EX','delta=e-4,t=0.1,EX','delta=e-4,t=0.5,EX',...
 %   'delta=e-3,t=0.1,UW','delta=e-3,t=0.5,UW','delta=e-4,t=0.1,UW','delta=e-4,t=0.5,UW')

% for i=1:20
%     plot(x0,data(i:i,3:13))
% end
%plot(x0,data(16:16,3:13),'r')

%plot(log(data(1:1,1:4)./data(1,4)),log(data(2:2,1:4)./data(2,4)),'r')
% plot(log(data(1:1,1:4)./data(1,4)),log(data(3:3,1:4)./data(3,4)),'b')
% plot(log(data(1:1,1:4)./data(1,4)),log(data(4:4,1:4)./data(4,4)),'g')
% plot(log(data(1:1,1:4)./data(1,4)),log(data(5:5,1:4)./data(5,4)),'k')
% plot(log(data(1:1,1:4)./data(1,4)),log(data(6:6,1:4)./data(6,4)),'y')
% plot(log(data(1:1,1:4)./data(1,4)),log(data(7:7,1:4)./data(7,4)),'m')
%legend('theta=0.5,mu=0.5')
% legend('theta=0,mu=0.3','theta=0,mu=0.5','theta=1,mu=0.3','theta=1,mu=0.5',...
%     'theta=0.5,mu=0.3','theta=0.5,mu=0.5')
