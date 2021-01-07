function Boost_DC(E,R,RL,L,f,D,C)
% programma BOOST_DC 06 gen 2021
% La function simula il comportamento della tensione di uscita v0 e della
% corrente nell'induttore iL pe run boost DC-DC
% La function richiama la function trattoC che definisce quel tratto
% notazioni:
% B - sw chiuso
% C - sw aperto
% v0 - ddp di uscita
% iL - corrente induttore; iLB corrente induttore nell'intervallo B...
%% PARAMETRI DI INPUT
% E amppiezza del segnale di ingresso [V]
% R resistenza di carico [Ohm]
% RL resistenza coil [Ohm]
% L induttore [H]
% f frequenza del segnale PWM [Hz]
% D duty cycle
% C capacità del condensatore [F]

VD=0.001; % tensione di soglia del diodo
T=1/f; % periodo
DT=D*T; % intervallo in cui l'interruttore è chiuso
t=linspace(0,T,1000); % vettore del tempo

% le funzioni sono diverse da zero solo nel loro intervallo
%% intervallo B (interruttore chiuso)
v0B=@(x,t) (t>=0).*(t<=DT).*(x.*exp(-1*(t)./(R.*C))); % ddp uscita

iLB=@(x,t) (t>=0).*(t<=DT).*((x-E/RL).*exp(-1*(t)./(L/RL))+ E/RL); %corrente induttore

%% imponiamo la continuità periodica
X=[(E-VD)/(R+RL),(E-VD)*R/(R+RL)]; % vettore delle condizioni iniziali I0 e V0
X0=[20 20];k=0;% vettori di appoggio
TT=linspace(DT,T,1000); % tempo di appoggio per il tratto C

while sum(abs(X-X0))>1e-5
    X0=X;k=k+1;
    
    I1=iLB(X(1),DT); % continuità
    V1=v0B(X(2),DT);
   
    [iLC,v0C]=trattoC([I1,V1],TT,E,R,RL,L,C,VD,T,DT);

    X(1)=iLC(end);
    X(2)=v0C(end);
    
    k,sum(abs(X-X0))
end % in uscita dal ciclo X contiene le c.i. periodiche

%% costruzione iL e v0
[iLC,v0C]=trattoC([I1,V1],t,E,R,RL,L,C,VD,T,DT);
iL=iLB(X(1),t)+iLC;
v0=v0B(X(2),t)+v0C;

%% PLOT
figure(1)
subplot(2,1,2),hold on, grid on, box on
plot(t+0.49,iL);grid on
legend('spice','matlab')
text(0,0,['L = ', num2str(L),' H ;']);


subplot(2,1,1),hold on, grid on, box on
plot(t+0.49,v0),grid on
legend('spice','matlab')

end
%% definizione funzioni intervallo C
function [iLC,v0C]=trattoC(X1,t,E,R,RL,L,C,VD,T,DT)
%global L C R RL E VD T DT

lam=roots([L*C L/R+RL*C RL/R+1]); % frequenze naturali

I1=X1(1); % continuità
V1=X1(2);
v0P=(E-VD)*R/(RL+R); % sol particolare costante
iLP=(E-VD)/(RL+R); % sol particolare costante

% calcolo costanti k1 e k2
A=[1 1;lam(1) lam(2)]; % matrice dei coefficienti
bL=@(I1,V1) [I1-iLP;(E-VD-V1-I1*RL)/L]; % termini noti per iL
k12=A\bL(I1,V1);
k1=k12(1);k2=k12(2);
%
bC=@(I1,V1) [V1-v0P;(I1-V1/R)/C]; % termini noti per v0
k34=A\bC(I1,V1);
k3=k34(1);k4=k34(2);

v0C=(t>=DT).*(t<=T).*...
    (k3.*exp(lam(1)*(t-DT))+k4.*exp(lam(2)*(t-DT))+v0P); % ddp uscita

iLC=(t>=DT).*(t<=T).*...
    (k1.*exp(lam(1)*(t-DT))+k2.*exp(lam(2)*(t-DT))+iLP); % induttore

I0=find(iLC<0, 1 ); % cerca il primo punto per cui il<0

if not(isempty(I0)) % se I0 non è vuoto (cioè ci sono punti in cui il<0)
    iLC(I0:end)=0*iLC(I0:end); %azzera iL
    v0C(I0:end)=v0C(I0)*exp(-1*(t(I0:end)-t(I0))./(R.*C)); %v0 evolve con costante di tempo RC
end
end
