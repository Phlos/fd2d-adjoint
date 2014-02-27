function f_signal=butterworth_lp(signal, t, npoles, cutoff, mode)

%- BUTTERWORTH LOWPASS FILTER ---------------------------------------------
%-
%- function f_signal=butterworth_lp(signal, t, npoles, cutoff, mode)
%-
%- signal = signal to be filtered
%- t = time vector
%- npoles = number of poles
%- cutoff = cutoff frequency in Hz
%- mode = 'silet' or 'plot'

n=npoles;
wc=2*pi*cutoff;

s=zeros(1,n);

t=t-t(1);

%- determine poles --------------------------------------------------------

i=sqrt(-1);

for k=1:n
    
    s(k)=wc*exp(i*pi*(2*k+n-1)/(2*n));
    
end

%- compute and plot transfer function -------------------------------------

x=-2*wc:(wc/100):2*wc;
y=-2*wc:(wc/100):2*wc;

[x y]=meshgrid(x,y);

H=ones(size(x));
P=ones(size(x));

for k=1:n
    
    P=P.*(x+i*y-s(k));
    
end

H=power(wc,n)./P;

if (strcmp(mode,'plot')==1)

    figure
    subplot(1,3,1);
    pcolor(x,y,abs(H));
    shading interp
    caxis([0 30])
    colormap(bone)
    phi=0:0.02:2.1*pi;
    xx=wc*cos(phi);
    yy=wc*sin(phi);
    hold on
    plot(xx,yy,'w');
    plot([0 0],[-2*wc 2*wc],'w');
    plot([-2*wc 2*wc],[0 0],'w');
    axis image

    subplot(1,3,2);
    plot(x(201,200:401)/(2*pi),abs(H(200:401,201)),'k');
    hold on
    plot([wc/(2*pi) wc/(2*pi)],[0 1.1],'b');
    axis([0 2*wc/(2*pi) 0 1.1]);

    subplot(1,3,3);
    plot(x(201,200:401)/(2*pi),unwrap(angle(H(200:401,201))),'k')
    hold on
    plot([wc/(2*pi) wc/(2*pi)],[-pi-0.1 0.1],'b');
    axis([0 2*wc/(2*pi) -pi-0.1 0.1]);

end
    
%- compute and plot impulse response --------------------------------------

Pi=ones(1,n);
h=0*t;

for i=1:n
    for k=1:n
        
        if (not(i==k))
            
            Pi(i)=Pi(i)*(s(i)-s(k));
            
        end
        
    end
end

Pi=Pi/power(wc,n);

for k=1:n
    
    h=h+exp(s(k)*t)/Pi(k);
    
end

if (strcmp(mode,'plot'))

    figure
    plot(t,real(h),'k');
    
end

f_signal=real(conv(h,signal))*(t(2)-t(1));
f_signal=f_signal(1:length(t));