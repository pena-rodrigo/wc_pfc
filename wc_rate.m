clear all

k=15; r0=0.5;
b=1;  x0=5;
% w=6.3; A = 2.35;
w=6; A=2.5;
r = -10:0.01:10;
planes=1;
Iinp=0.8; %osc input amp
transient = 1000; %ms

tf=2000+transient; dt=0.1; N=(tf/dt);
ttstart = 0.1;
ttstop  = tf;
fi=1; ff=50; % [Hz]
df=1;

taur=1; taun=100;

Rinf1 = @(x) 1./(1+exp(-(x-x0)));
Rinf2 = @(x) 1./(1+exp(-(x-1.2*x0)));
ainf = @(r) 1./(1+exp(-k*(r-r0)));
binf = @(r) 1./(1+exp(k*(r-r0)));
cinf = @(r) 1./(1+exp(-k*(r-r0)));
dinf = @(r) 1./(1+exp(k*(r-r0)));
anull = @(r,A) (-w*r - A + x0 + log(r./(1-r)))/b;

if(planes==1)
    figure(1);
    b=-1;
    anull = @(r,A) (-w*r - A + x0 + log(r./(1-r)))/b;
    subplot(2,2,1)
    plot(r,ainf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
%     plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
%     plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    set(gca,'xtick',[]);
    ylabel('second-variable')
    title('negative feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    subplot(2,2,2)
    plot(r,binf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
%     plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
%     plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    set(gca,'xtick',[]); set(gca,'ytick',[])
    title('positive feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    b=1;
    x0=1.2*x0;
    anull = @(r,A) (-w*r - A + x0 + log(r./(1-r)))/b; %b=-1 and x0 = 1.2
    Rinf = @(x) 1./(1+exp(-(x-x0)));
    subplot(2,2,3)
    plot(r,cinf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
%     plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
%     plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    xlabel('r')
    ylabel('second-variable')
    title('positive feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

    subplot(2,2,4)
    plot(r,dinf(r),'g','LineWidth',2); hold on;
    plot(r,anull(r,A),'r','LineWidth',2); hold on;
%     plot(r,anull(r,A+Iinp),'r--','LineWidth',2); hold on;
%     plot(r,anull(r,A-Iinp),'r--','LineWidth',2); hold on;
    xlim([0,1]); ylim([0,1]);
    xlabel('r')
    set(gca,'ytick',[])
    title('negative feedback')
    set(gca, 'FontName', 'Times New Roman','FontSize',20)

end

b=1; w=6; A = 2.5; x0=5;
acum1=[];
acum2=acum1; acum3=acum1; acum4=acum1;
for f = fi:df:ff
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
    Input=zeros(1,N);

    for i = 1:N-1
        t = i*dt;
        if (t<transient)
            I=A;
        else
            I = A+Iinp*sin(2*pi*f*t*10^-3);
        end
       
        %a
        r1(i+1) = r1(i) + dt*(-r1(i) + Rinf1(w*r1(i) - b*n1(i) + I +r2(i)))/taur;
        n1(i+1) = n1(i) + dt*(-n1(i) + ainf(r1(i)))/taun;
        %b
        r2(i+1) = r2(i) + dt*(-r2(i) + Rinf1(w*r2(i) - b*n2(i) + I))/taur;
        n2(i+1) = n2(i) + dt*(-n2(i) + binf(r2(i)))/taun;
        %c
        r3(i+1) = r3(i) + dt*(-r3(i) + Rinf2(w*r3(i) + b*n3(i) + I ))/taur;
        n3(i+1) = n3(i) + dt*(-n3(i) + cinf(r3(i)))/taun;
        %d
        r4(i+1) = r4(i) + dt*(-r4(i) + Rinf2(w*r4(i) + b*n4(i) + I +r3(i)))/taur;
        n4(i+1) = n4(i) + dt*(-n4(i) + dinf(r4(i)))/taun;          
        Input(i) = I;
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
    
    if(f==1)
        figure(1)
        subplot(2,2,1)
        plot(r1(1:end),n1(1:end),'-','Color',[0,1,1],'LineWidth',0.1) %end/2
        subplot(2,2,2)
        plot(r2(1:end),n2(1:end),'-','Color',[0,1,1],'LineWidth',0.1)
        subplot(2,2,3)
        plot(r3(1:end),n3(1:end),'-','Color',[0,1,1],'LineWidth',0.1)
        subplot(2,2,4)
        plot(r4(1:end),n4(1:end),'-','Color',[0,1,1],'LineWidth',0.1)  
        
        figure(4);
        subplot(2,2,1)
        plot(dt:dt:dt*N,r1,'Color',[0,1,1],'LineWidth',1); hold on;
        set(gca,'xtick',[]);
        ylabel('r')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,2)
        plot(dt:dt:dt*N,r2,'Color',[0,1,1],'LineWidth',1); hold on;
        set(gca,'xtick',[]); set(gca,'ytick',[])
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,3)
        plot(dt:dt:dt*N,r3,'Color',[0,1,1],'LineWidth',1); hold on;
        ylabel('r')
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,4)
        plot(dt:dt:dt*N,r4,'Color',[0,1,1],'LineWidth',2); hold on;
        set(gca,'ytick',[])
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
    elseif(f==4)
        figure(1)
        subplot(2,2,1)
        plot(r1(1:end),n1(1:end),'-','Color',[1,0,1],'LineWidth',0.1)
        subplot(2,2,2)
        plot(r2(1:end),n2(1:end),'-','Color',[1,0,1],'LineWidth',0.1)
        subplot(2,2,3)
        plot(r3(1:end),n3(1:end),'-','Color',[1,0,1],'LineWidth',0.1)
        subplot(2,2,4)
        plot(r4(1:end),n4(1:end),'-','Color',[1,0,1],'LineWidth',0.1)  
        
        figure(4);
        subplot(2,2,1)
        plot(dt:dt:dt*N,r1,'Color',[1,0,1],'LineWidth',1)
        set(gca,'xtick',[]);
        ylabel('r')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,2)
        plot(dt:dt:dt*N,r2,'Color',[1,0,1],'LineWidth',1)
        set(gca,'xtick',[]); set(gca,'ytick',[])
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,3)
        plot(dt:dt:dt*N,r3,'Color',[1,0,1],'LineWidth',1)
        ylabel('r')
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,4)
        plot(dt:dt:dt*N,r4,'Color',[1,0,1],'LineWidth',1)
        set(gca,'ytick',[])
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
    elseif(f==50)
        figure(1)
        subplot(2,2,1)
        plot(r1(1:end),n1(1:end),'-','Color',[0,0,1],'LineWidth',0.1)
        subplot(2,2,2)
        plot(r2(1:end),n2(1:end),'-','Color',[0,0,1],'LineWidth',0.1)
        subplot(2,2,3)
        plot(r3(1:end),n3(1:end),'-','Color',[0,0,1],'LineWidth',0.1)
        subplot(2,2,4)
        plot(r4(1:end),n4(1:end),'-','Color',[0,0,1],'LineWidth',0.1)  
        
        figure(4);
        subplot(2,2,1)
        plot(dt:dt:dt*N,r1,'Color',[0,0,1],'LineWidth',1)
        set(gca,'xtick',[]);
        ylabel('r')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,2)
        plot(dt:dt:dt*N,r2,'Color',[0,0,1],'LineWidth',1)
        set(gca,'xtick',[]); set(gca,'ytick',[])
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,3)
        plot(dt:dt:dt*N,r3,'Color',[0,0,1],'LineWidth',1)
        ylabel('r')
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        subplot(2,2,4)
        plot(dt:dt:dt*N,r4,'Color',[0,0,1],'LineWidth',1)
        set(gca,'ytick',[])
        xlabel('Time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
    end
    
    if(f==1)
        figure(10)
        subplot(1,3,1)
        plot(dt:dt:tf-dt,Input(1:end-1),'Color',[0,1,1],'LineWidth',2); 
        title('frequency = 1Hz')
        ylabel('input')
        xlabel('time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        xlim([0,tf])
    elseif(f==10)
        figure(10)
        subplot(1,3,2)
        plot(dt:dt:tf-dt,Input(1:end-1),'Color',[1,0,1],'LineWidth',2)
        title('frequency = 10Hz')
        xlabel('time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)    
        xlim([0,tf])
    elseif(f==50)
        figure(10)
        subplot(1,3,3)
        plot(dt:dt:tf-dt,Input(1:end-1),'Color',[0,0,1],'LineWidth',2)
        title('frequency = 50Hz')
        xlabel('time [ms]')
        set(gca, 'FontName', 'Times New Roman','FontSize',20)
        xlim([0,tf])
    end
    
    pks1 = [pks1; length(findpeaks(r1))/2];
    pks2 = [pks2; length(findpeaks(r4))/2];
    
    acum1 = [acum1; max(r1(end/2:end)) - min(r1(end/2:end))];% sum(r1)];
    acum2 = [acum2; max(r2(end/2:end)) - min(r2(end/2:end))];% sum(r2)];
    acum3 = [acum3; max(r3(end/2:end)) - min(r3(end/2:end))];% sum(r3)];
    acum4 = [acum4; max(r4(end/2:end)) - min(r4(end/2:end))];% sum(r4)];
end



figure; 
subplot(2,1,1)
plot(pks1)
xlabel('\tau_a [ms]')
ylabel('Frequency [Hz]')
title('1')
set(gca, 'FontName', 'Times New Roman','FontSize',14)
subplot(2,1,2)
plot(pks2)
xlabel('\tau_a [ms]')
ylabel('Frequency [Hz]')
title('4')
set(gca, 'FontName', 'Times New Roman','FontSize',14)

figure(3);
subplot(2,2,1)
plot(acum1,'k','LineWidth',2); hold on;
plot(1,acum1(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum1(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum1(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('r_{\rm max} - r_{\rm min}')
set(gca, 'FontName', 'Times New Roman','FontSize',14)

subplot(2,2,2)
plot(acum2,'k','LineWidth',2); hold on;
plot(1,acum2(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum2(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum2(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
set(gca, 'FontName', 'Times New Roman','FontSize',14)

subplot(2,2,3)
plot(acum3,'k','LineWidth',2); hold on;
plot(1,acum3(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum3(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum3(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
ylabel('r_{\rm max} - r_{\rm min}')
xlabel('frequency [Hz]')
set(gca, 'FontName', 'Times New Roman','FontSize',14)

subplot(2,2,4)
plot(acum4,'k','LineWidth',2); hold on;
plot(1,acum4(1),'o','Color',[0,1,1],'MarkerFaceColor',[0,1,1])
plot(10,acum4(10),'o','Color',[1,0,1],'MarkerFaceColor',[1,0,1])
plot(50,acum4(50),'o','Color',[0,0,1],'MarkerFaceColor',[0,0,1])
xlabel('frequency [Hz]')
set(gca, 'FontName', 'Times New Roman','FontSize',14)