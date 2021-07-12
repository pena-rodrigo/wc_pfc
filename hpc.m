clear all
% close all

k=15; r0=0.5;
b=2;  x0=5;
% w=6.3; A = 2.35;
w=6; A=2.5; 
r = -10:0.01:10;
planes=0;
Iinp=1;0.8; %osc input amp
transient = 1000; %ms

tf=2000+transient; dt=0.1; N=(tf/dt);
ttstart = 0.1;
ttstop  = tf;
fi=1; ff=50; % [Hz]
df=1;

D=500e-4;  
count=1;
for D = 0:10:50
D = D*1e-4

taur=1; taun=25;

Rinf1 = @(x) 1./(1+exp(-(x-x0)));
ainf = @(r) 1./(1+exp(-k*(r-r0)));
anull = @(r,A) -(-w*r - A + x0 + log(r./(1-r)))/b;

if(planes==1)
    figure(1);
    subplot(2,2,1)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    ylabel('second-variable')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    subplot(2,2,2)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    set(gca, 'FontName', 'Times New Roman','FontSize',20)
    
    subplot(2,2,3)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    
    subplot(2,2,4)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    set(gca, 'FontName', 'Times New Roman','FontSize',20)


end

% b=1; w=6; A = 2.5; x0=5;
acum1=[];
acum2=acum1; acum3=acum1; acum4=acum1; acum5=acum1; acum6=acum1;
f2 = 15;
for f1 = 1:50 
% f2 = f1;

% f=0; 
pks1=[]; pks2=[];
% for taun = 100:1:100
    r1=zeros(1,N); r1(1)=0.4;
    n1=zeros(1,N); n1(1)=0.5;
    r2=zeros(1,N); r2(1)=0.4;
    n2=zeros(1,N); n2(1)=0.5;
    r3=zeros(1,N); r3(1)=0.4;
    n3=zeros(1,N); n3(1)=0.5;
    r4=zeros(1,N); r4(1)=0.4;
    n4=zeros(1,N); n4(1)=0.5;
    r5=zeros(1,N); r5(1)=0.4;
    n5=zeros(1,N); n5(1)=0.5;
    r6=zeros(1,N); r6(1)=0.4;
    
    Input=zeros(1,N);

    for i = 1:N-1
        t = i*dt;
        if (t<transient)
            I1=A;
            I2=A;           
        else
            I1 = A+Iinp*sin(2*pi*f1*t*10^-3);
            I2 = A+Iinp*sin(2*pi*f2*t*10^-3);
        end
       
        r1(i+1) = r1(i) + dt*(-r1(i) + Rinf1(w*r1(i) - b*n1(i) + I1))/taur;
        r1(i+1) = r1(i+1) + sqrt(dt*2*D)*rand()/taur;
        
        n1(i+1) = n1(i) + dt*(-n1(i) + ainf(r1(i)+r2(i)))/taun;
        n1(i+1) = n1(i+1) + sqrt(dt*2*D)*rand()/taun;
        
        r2(i+1) = r2(i) + dt*(-r2(i) + Rinf1(w*r2(i) - b*n1(i) - b*n2(i) + I1))/taur;    
        r2(i+1) = r2(i+1) + sqrt(dt*2*D)*rand()/taur;
        
        n2(i+1) = n2(i) + dt*(-n2(i) + ainf(r2(i)+r3(i)))/taun;
        n2(i+1) = n2(i+1) + sqrt(dt*2*D)*rand()/taun;
        
        r3(i+1) = r3(i) + dt*(-r3(i) + Rinf1(w*r3(i) - b*n2(i) - b*n3(i) + I1))/taur;
        r3(i+1) = r3(i+1) + sqrt(dt*2*D)*rand()/taur;
        
        n3(i+1) = n3(i) + dt*(-n3(i) + ainf(r3(i)+r4(i)))/taun;
        n3(i+1) = n3(i+1) + sqrt(dt*2*D)*rand()/taun;
        
        r4(i+1) = r4(i) + dt*(-r4(i) + Rinf1(w*r4(i) - b*n3(i) - b*n4(i) + I2))/taur; 
        r4(i+1) = r4(i+1) + sqrt(dt*2*D)*rand()/taur;
        
        n4(i+1) = n4(i) + dt*(-n4(i) + ainf(r4(i)+r5(i)))/taun;
        n4(i+1) = n4(i+1) + sqrt(dt*2*D)*rand()/taun;
        
        r5(i+1) = r5(i) + dt*(-r5(i) + Rinf1(w*r5(i) - b*n4(i) - b*n5(i) + I2))/taur; 
        r5(i+1) = r5(i+1) + sqrt(dt*2*D)*rand()/taur;
        
        n5(i+1) = n5(i) + dt*(-n5(i) + ainf(r5(i)+r6(i)))/taun;
        n5(i+1) = n5(i+1) + sqrt(dt*2*D)*rand()/taun;
        
        r6(i+1) = r6(i) + dt*(-r6(i) + Rinf1(w*r6(i) - b*n5(i) + I2))/taur; 
        r6(i+1) = r6(i+1) + sqrt(dt*2*D)*rand()/taur;
    
%         Input(i) = I;
    end
    
%     figure(1)
%     subplot(2,2,1)
%     plot(r1,n1,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,2)
%     plot(r2,n2/2+n1/2,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,3)
%     plot(r3,n2/2+n3/2,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,4)
%     plot(r4,n3,'-','Color',[0,1,1],'LineWidth',0.1)
    
%         figure(1)
%         subplot(1,2,1)
%         plot(r1(1:end),n1(1:end),'-','Color',[0,1,1],'LineWidth',0.1) %end/2
%         subplot(1,2,2)
%         plot(r2(1:end),n1(1:end),'-','Color',[0,1,1],'LineWidth',0.1)
% 
%         figure(4);
%         subplot(1,2,1)
%         plot(dt:dt:dt*N,r1,'Color',[0,1,1],'LineWidth',1); hold on;
%         set(gca,'xtick',[]);
%         ylabel('r')
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)
%         subplot(1,2,2)
%         plot(dt:dt:dt*N,r2,'Color',[0,1,1],'LineWidth',1); hold on;
%         set(gca,'xtick',[]); set(gca,'ytick',[])
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)
%       
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)
%   
    
%     if(f==1)
%         figure(10)
%         subplot(1,3,1)
%         plot(dt:dt:tf-dt,Input(1:end-1),'Color',[0,1,1],'LineWidth',2); 
%         title('frequency = 1Hz')
%         ylabel('input')
%         xlabel('time [ms]')
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)
%         xlim([0,tf])
%     elseif(f==10)
%         figure(10)
%         subplot(1,3,2)
%         plot(dt:dt:tf-dt,Input(1:end-1),'Color',[1,0,1],'LineWidth',2)
%         title('frequency = 10Hz')
%         xlabel('time [ms]')
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)    
%         xlim([0,tf])
%     elseif(f==50)
%         figure(10)
%         subplot(1,3,3)
%         plot(dt:dt:tf-dt,Input(1:end-1),'Color',[0,0,1],'LineWidth',2)
%         title('frequency = 50Hz')
%         xlabel('time [ms]')
%         set(gca, 'FontName', 'Times New Roman','FontSize',20)
%         xlim([0,tf])
%     end
    
    pks1 = [pks1; length(findpeaks(r1))/2];
    pks2 = [pks2; length(findpeaks(r2))/2];
    
    acum1 = [acum1; max(r1(end/2:end)) - min(r1(end/2:end))];% sum(r1)];
    acum2 = [acum2; max(r2(end/2:end)) - min(r2(end/2:end))];% sum(r2)];
    acum3 = [acum3; max(r3(end/2:end)) - min(r3(end/2:end))];% sum(r1)];
    acum4 = [acum4; max(r4(end/2:end)) - min(r4(end/2:end))];% sum(r2)];
    acum5 = [acum5; max(r5(end/2:end)) - min(r5(end/2:end))];% sum(r1)];
    acum6 = [acum6; max(r6(end/2:end)) - min(r6(end/2:end))];% sum(r2)];
    
end
% figure; plot(sum(r1)); hold on; plot(sum(r2))
% figure; 
% subplot(2,1,1)
% plot(pks1)
% xlabel('\tau_a [ms]')
% ylabel('Frequency [Hz]')
% title('1')
% set(gca, 'FontName', 'Times New Roman','FontSize',14)
% subplot(2,1,2)
% plot(pks2)
% xlabel('\tau_a [ms]')
% ylabel('Frequency [Hz]')
% title('4')
% set(gca, 'FontName', 'Times New Roman','FontSize',14)

% figure(3);
% subplot(1,2,1)
% plot(acum1,'k','LineWidth',2); hold on;
% plot(1,acum1(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
% plot(10,acum1(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
% plot(50,acum1(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
% ylabel('r_{\rm max} - r_{\rm min}')
% set(gca, 'FontName', 'Times New Roman','FontSize',14)
% 
% subplot(1,2,2)
% plot(acum2,'k','LineWidth',2); hold on;
% plot(1,acum2(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
% plot(10,acum2(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
% plot(50,acum2(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
% set(gca, 'FontName', 'Times New Roman','FontSize',14)
% 
% 
figure; 
plot(smooth(acum1),'b','LineWidth',2); 
hold on; 
plot(smooth(acum2),'r','LineWidth',2);
plot(smooth(acum3),'g','LineWidth',2); 
plot(smooth(acum4),'k','LineWidth',2);
plot(smooth(acum5),'c','LineWidth',2); 
plot(smooth(acum6),'color',[1,0,1],'LineWidth',2);
legend('1','2','3','4','5','6')

figure(); 
subplot(6,1,1);
plot(dt:dt:tf,r1,'LineWidth',1);
ylabel('r_1')
set(gca,'xtick',[]);
set(gca, 'FontName', 'Helvetica','FontSize',20)

subplot(6,1,2);
plot(dt:dt:tf,r2,'LineWidth',1);
ylabel('r_2')
set(gca,'xtick',[]);
set(gca, 'FontName', 'Helvetica','FontSize',20)

subplot(6,1,3);
plot(dt:dt:tf,r3,'LineWidth',1);
ylabel('r_3')
set(gca,'xtick',[]);
set(gca, 'FontName', 'Helvetica','FontSize',20)

subplot(6,1,4);
plot(dt:dt:tf,r4,'LineWidth',1); 
ylabel('r_4')
set(gca,'xtick',[]);
set(gca, 'FontName', 'Helvetica','FontSize',20)

subplot(6,1,5); 
plot(dt:dt:tf,r5,'LineWidth',1); 
set(gca,'xtick',[]);
ylabel('r_5')
set(gca, 'FontName', 'Helvetica','FontSize',20)

subplot(6,1,6);
plot(dt:dt:tf,r6,'LineWidth',1);
xlabel('Time (ms)');
ylabel('r_6')
set(gca, 'FontName', 'Helvetica','FontSize',20)

figure(15); 
subplot(2,3,count)
hold on;
plot(smooth((acum1+acum2+acum3)-(acum4+acum5+acum6)),'LineWidth',2,'Color',[0 0 .5]);
plot(1:50,zeros(1,50),'--','LineWidth',1,'Color',[.3 .3 .3]) %[.75 .75 1]
legend(['D=' num2str(D)]);
xlim([1,50])
ylim([-0.8,0.65])
if(count>3)
    xlabel('Input frequency (Hz)')
end
ylabel('Selection index')
set(gca, 'FontName', 'Helvetica','FontSize',20)
count=count+1;

% 
end