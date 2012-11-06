function [x,y] = diskrettime( NT,N,T )
%Loesung der Waermeleitungsgleichung in der Zeit,
%diskret in Ort und Raum.
%NT - Zeitschritte pro Sekunde, N - Raumschritte,
%T - Endzeit. Achtung: Stabilitaetsbedingung beachten!

if (nargin<3)
T=1;
end
if (nargin<1)
NT=1000;
N=70; %Ab einem Wert grüßer als 70 wird Chaos erzeugt
end
NT=round(T*NT);
disp=round(NT/100);
dt=T/NT;
h=pi/N;
%Stabilitätskonstante (<1/2)
%dt/(h*h)
x=((1:N-1)*h)';
A=(-2*diag(ones(N-1,1))+diag(ones(N-2,1),1)+diag(ones(N-2,1),-1))/(h*h);
y=zeros(size(x));
for i=1:NT
    y=y+dt*(A*y+q(x));
    if (mod(i,disp)==1)
    plot(x,y);
    t=i*dt;
    title (['Temperaturverlauf zum Zeitpunkt t=' num2str(t)]);
    drawnow;    
    end
end
end
function g=q(x)
g=zeros(size(x));
g(abs(x-pi/2)<0.5)=1;
end
