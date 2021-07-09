%definir flujo volumétrico
Q=0.00726496; %m3/s;
a1=9.97e-5;%1/s
a2=1.07e-3;%1/s
a3=5.70e5;%1/s
g1=17254.4;%K
g2=10224.4;%K
l1=-9061227;%K^2
l2=-6123081;%K^2
e3=8.1e1;%kJ/mol
delta=1.45;
Tk=623;%K
R=8.314e-3;%kJ/mol K
k1=a1*exp((g1/Tk)+(l1/(Tk^2))); %1/s
k2=a2*exp((g2/Tk)+(l2/(Tk^2))); %1/s
k3=a3*exp(-e3/(R*Tk)); %1/s
B0=1;%kgi/kgT
C0=0;%kgi/kgT
W0=1;%kgi/kgT
t0=0;%s
h=1e-3;%adimensional
B(1)=B0;%kgi/kgT
C1(1)=C0;%kgi/kgT
C2(1)=0;%kgi/kgT
G(1)=0;%kgi/kgT
W(1)=W0;%kgi/kgT
t(1)=t0; %s
p=650;%kgTs/m3
e=1-0.1;%adimensional
tao(1)=0;%s
Weje=1571;%W/m3
Hrx=-235;%kJ/kg %%%teniendo en cuenta el efecto de la temperatura y el tamaño de partícula sobre la pirólisis (Koufopanos et al. 1991)
cal(1)=0;%kJ/kg
for i=1:600000
t(i+1)=t(i)+h;
dB(i)=-1*(k1+k2)*B(i);
B(i+1)=B(i)+h*dB(i);
dC1(i)=k2*B(i)-k3*C1(i);
C1(i+1)=C1(i)+h*dC1(i);
dC2(i)=delta*k3*C1(i);
C2(i+1)=C2(i)+h*dC2(i);
dW(i)=-1*k1*B(i)-k3*C1(i)+delta*k3*C1(i);
W(i+1)=W(i)+h*dW(i);
dG(i)=k1*B(i);
G(i+1)=G(i)+h*dG(i);
dGc(i)=k1*B(i)-delta*k3*C1(i)-k2*B(i)-k3*C1(i);
dtotal(i)=-k1*B(i)-k2*B(i)-k3*C1(i)*G(i);%rtotal
%Balances de materia
%Fase sólida
%entrada-salida+generación=acumulación
%entrada=Qs(i)Cixs
%salida=Qs(i)Cix+xds
%generación=dV*ri=Alecho*dx*ri
%tras hacer las consideraciones y operar los factores, se obtiene:
%dcixs=p*(1-e)*ri(s)
%i corresponde a cada uno de los productos y reactivos de la pirólisis
tao(i+1)=tao(i)+h;
Bx(i+1)=p*(1-e)*B(i);%kg/m3
C1x(i+1)= p*(1-e)*C1(i);%kg/m3
C2x(i+1)= p*(1-e)*C2(i);%kg/m3
Wx(i+1)=p*(1-e)*W(i);%kg/m3
%fase gaseosa 
Gx(i+1)=p*(1-e)*G(i);%kg/m3
%Balances de energía
dq(i+1)= -1*(-p*(1-e)*dtotal(i)*Hrx);%kW/m3
%cal(i+1)=cal(i)+h*dq(i);
end 
%flujos másicos de salida%
%relación de concentraciones
rel=.10929;
biomasa=Bx(600001)*(Q/rel);%kg/s
char1=C1x(600001)*(Q/rel);%kg/s
char2=C2x(600001)*(Q/rel);%kg/s
gasesyvolatiles=Gx(600001)*(Q/rel);%kg/s
%flujo de nitrógeno como atmósfera inerte
vnitroteor=2.5;%cm3/min por cada 20mg de biomasa. [koufopanos 1989]
vnitro=((vnitroteor*(1/60)*(1/1e+6)*(Q*p)*1000)/(20e-3));%m3/s
pnitro=1.126;%kg/m3,https://www.engineeringtoolbox.cofv.m/nitrogen-N2-density-specific-weight-temperature-pressure-d_2039.html 
mnitro=vnitro*pnitro;
%%%fijar el diametro%%%
D=1;%m
%%%velocidad de rotación del horno%%%
omega=2.00;%rpm
%%%tiempo máximo de residencia%%%
tau=600;%s
%%%volumen del reactor%%%
V=Q*tau;%m3
%%%área de sección transversal del reactor%%%
A=pi()*(D^2)/4;%m2
%%%longitud del reactor%%%
L=V/A;%m
%%%graficos generados, para ejecutarlos por favor borre el símbolo de porcentaje que se presenta antes de cada plot%%%
%plot(tao,C1x)
%plot(tao,C2x)
plot(tao,Wx)
hold on
plot(tao,Gx)
%axis([0 650 0 70]);
%plot(tao,dq)
%plot(tao,cal)
qpiro=trapz(tao,dq);%demanda energética por unidad de volúmen kW/m3
qtotal=qpiro*V;%demanda energética kW
%%%Fluido de calentamiento diesel%%%
%ecuación de Aly-lee NIST para el cálculo de cp recuperado de bases de
%datos Aspen
c1i=248976.8;%J/kmolK
c2i=77686.8;%J/kmolK
c3i=1235.959;%K
c4i=400927.5;%J/kmolK
c5i=614.3928;%K
%cpig=c1i+c2i*((c3i/Temp)/sinh(c3i/Temp))^2+c4i*((c5i/Temp)/cosh(c5i/Temp))^2
syms x    
fun= c1i+(c2i*((c3i/x)/sinh(c3i/x))^2)+(c4i*((c5i/x)/cosh(c5i/x))^2);%fuentes Aspen Plus V11.
q=int(fun,640,800); %fuentes Aspen Plus V11.
hd2=double(q)*(1/1000)*(1/1000);%kJ/mol %fuentes Aspen Plus V11.
md2=((qtotal)/(hd2))*(226.41)*(1/1000);%kg/s %fuentes Aspen Plus V11.
vd2=md2/344.44;%m3/s %fuentes Aspen Plus V11.
%evaluado a 573 623 673 873