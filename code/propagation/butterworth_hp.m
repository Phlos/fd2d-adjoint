function f_signal=butterworth_hp(signal, t, npoles, cutoff, mode)

%- BUTTERWORTH HIGHPASS FILTER --------------------------------------------
%-
%- function f_signal=butterworth_lp(signal, t, npoles, cutoff, mode)
%-
%- signal = signal to be filtered
%- t = time vector
%- npoles = number of poles
%- cutoff = cutoff frequency in Hz
%- mode = 'silent' or 'plot'

n=npoles;
wc=2*pi*cutoff;

s=zeros(1,n);

t=t-t(1);

%- determine poles --------------------------------------------------------

i=sqrt(-1);

for k=1:n
    
    s(k)=-wc*exp(i*pi*(n+1-2*k)/(2*n));
    
end

%- compute and plot transfer function -------------------------------------

x=-4*wc:(wc/100):4*wc;
y=-4*wc:(wc/100):4*wc;

[x y]=meshgrid(x,y);

H=ones(size(x));
P=ones(size(x));

for k=1:n
    
    P=P.*(1./(x+i*y+eps)-1/s(k));
    
end

H=1./(power(wc,n)*P);

if (strcmp(mode,'plot')==1)

    figure
    subplot(1,3,1);
    pcolor(x,y,abs(H));
    shading interp
    caxis([0 15])
    colormap(bone)
    phi=0:0.02:2.1*pi;
    xx=wc*cos(phi);
    yy=wc*sin(phi);
    hold on
    plot(xx,yy,'w');
    plot([0 0],[-2*wc 2*wc],'w');
    plot([-2*wc 2*wc],[0 0],'w');
    axis image

    size(x)
    
    subplot(1,3,2);
    plot(x(401,:)/(2*pi),abs(H(:,401)),'k');
    hold on
    plot([wc/(2*pi) wc/(2*pi)],[-0.1 1.1],'b');
    %axis([0 2*wc/(2*pi) 0 1.1]);

    subplot(1,3,3);
    plot(x(401,:)/(2*pi),(angle(H(:,401))),'k')
    hold on
    plot([wc/(2*pi) wc/(2*pi)],[-pi-0.1 pi+0.1],'b');
    %axis([0 2*wc/(2*pi) -pi-0.1 pi+0.1]);

end
    


%- compute and plot impulse response --------------------------------------

Pi=ones(1,n);
h=0*t;
g=0*t;

for i=1:n
    for k=1:n
        
        if (not(i==k))
            
            Pi(i)=Pi(i)*(1-s(i)/s(k));
            
        end
        
    end
end

Pi=Pi*power(wc,n);

for k=1:n
    
    h=h-power(s(k),n)*exp(s(k)*t)/Pi(k);
    
end

%h(1)=h(1)/(t(2)-t(1));

for k=2:length(t)-1
    h(k)=(h(k+1)-h(k))/(t(2)-t(1));
end

h(1)=h(1)-sum(h);
m=max(abs(fft(h)))*(t(2)-t(1));
h=h/m;

if (strcmp(mode,'plot'))

    figure
    plot(t,real(h),'k');
    
end

f_signal=real(conv(h,signal))*(t(2)-t(1));
f_signal=f_signal(1:length(t));
