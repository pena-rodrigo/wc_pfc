clear all
close all

k=15; r0=0.5;
b=1;  x0=5;
% w=6.3; A = 2.35;
w=6; A=2.5; 
r = -10:0.01:10;
planes=1;
Iinp=1;0.8; %osc input amp
transient = 1000; %ms

tf=2000+transient; dt=0.1; N=(tf/dt);
ttstart = 0.1;
ttstop  = tf;
fi=1; ff=50; % [Hz]
df=1;
% D=0.001;
rateP=9000;
weight=0.1;

taur=1; taun=25;

Rinf1 = @(x) 1./(1+exp(-(x-x0)));
ainf = @(r) 1./(1+exp(-k*(r-r0)));
anull = @(r,A) -(-w*r - A + x0 + log(r./(1-r)))/b;

if(planes==1)
    figure(1);
    subplot(1,2,1)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    ylabel('second-variable')
    title('negative feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    subplot(1,2,2)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
    plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
    plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    title('positive feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

end

% b=1; w=6; A = 2.5; x0=5;
acum1=[];
acum2=acum1; acum3=acum1; acum4=acum1;
f2 = 20;
for f1 = 1:50
% f2 = f1;

% f=0; 
pks1=[]; pks2=[];
% for taun = 100:1:100
    r1=zeros(1,N); r1(1)=0.4;
    n1=zeros(1,N); n1(1)=0.5;
    r2=zeros(1,N); r2(1)=0.4;
    n2=zeros(1,N); n2(1)=0.5;

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
%         r1(i+1) = r1(i+1) + sqrt(dt*2*D)*rand()/taur;
        if(rand()<(rateP*dt)/1000) %divide by 1000 to adjust units [kHz]
            r1(i+1) = r1(i+1) + weight*dt/taur;
        end
        
        n1(i+1) = n1(i) + dt*(-n1(i) + ainf(r1(i)+r2(i)))/taun;
       
        r2(i+1) = r2(i) + dt*(-r2(i) + Rinf1(w*r2(i) - b*n1(i) + I2))/taur;
%         r2(i+1) = r2(i+1) + sqrt(dt*2*D)*rand()/taur;
        if(rand()<(rateP*dt)/1000) %divide by 1000 to adjust units [kHz]
            r2(i+1) = r2(i+1) + weight*dt/taur;
        end
        

    
%         Input(i) = I;
    end
    
    figure(1)
%     subplot(2,2,1)
%     plot(r1,n1,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,2)
%     plot(r2,n2,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,3)
%     plot(r3,n3,'-','Color',[0,1,1],'LineWidth',0.1)
%     subplot(2,2,4)
%     plot(r4,n4,'-','Color',[0,1,1],'LineWidth',0.1)
    
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

figure(3);
subplot(1,2,1)
plot(acum1,'k','LineWidth',2); hold on;
plot(1,acum1(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum1(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum1(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('r_{\rm max} - r_{\rm min}')
set(gca, 'FontName', 'Times New Roman','FontSize',14)

subplot(1,2,2)
plot(acum2,'k','LineWidth',2); hold on;
plot(1,acum2(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum2(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum2(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
set(gca, 'FontName', 'Times New Roman','FontSize',14)


figure; plot(acum1,'b','LineWidth',2); 
hold on; 
plot(acum2,'r','LineWidth',2);
legend('variable','f=20Hz')


% 