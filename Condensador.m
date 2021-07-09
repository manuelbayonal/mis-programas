%1 acido ac�tico
%2 butanal
%3 acido but�rico
%4 alcohol conifer�lico
%5 �cido f�rmico
%6 guaiacol
%7 pentanal 
%8 fenol
%9 propanal
%10 �cido propanoico
%11 vapor de agua
%12 nitr�geno
F=(gasesyvolatiles/(Gx(600001)*10))+vnitro;%Flujo volum�trico
p1=351;%kg/m3
p2=419.07;%kg/m3
p3=469.41;%kg/m3
p4=839.09;%kg/m3
p5=1207.3;%kg/m3
p6=732.1;%kg/m3
p7=397.98;%kg/m3
p8=399.43;%kg/m3
p9=703.52;%kg/m3
p10=487.82;%kg/m3
p11=580.99;%kg/m3
p12=1.126;%kg/m3
x1=0.0370/(1+0.0288);
x2=0.109/(1+0.0288);
x3=0.011/(1+0.0288);
x4=0.190/(1+0.0288);
x5=0.042/(1+0.0288);
x6=0.108/(1+0.0288);
x7=0.021/(1+0.0288);
x8=0.054/(1+0.0288);
x9=0.144/(1+0.0288);
x10=0.017/(1+0.0288);
x11=0.267/(1+0.0288);
x12=0.0288/(1+0.0288);
R=8.314; %J/mol*K
pc1=57.10; %atm
pc2=40.00; %atm
pc3=52.00; %atm
pc4=33.60; %atm
pc5=57.34; %atm
pc6=46.61; %atm
pc7=35.00; %atm
pc8=60.50; %atm
pc9=47.00; %atm
pc10=53.00; %atm
pc11=217.60; %atm
pc12=33.53; %atm
tc1=594.00; %K
tc2=524.00; %K
tc3=628.00; %K
tc4=569.90; %K
tc5=580.00; %K
tc6=696.80; %K
tc7=554.00; %K
tc8=694.20; %K
tc9=496.00; %K
tc10=612.00; %K
tc11=647.20; %K
tc12=126.19; %K
tb1=391.15; %K
tb2=347.95; %K
tb3=436.65; %K
tb4=437.15; %K
tb5=373.95; %K
tb6=478.15; %K
tb7=411.15; %K
tb8=454.85; %K
tb9=321.95; %K
tb10=414.35; %K
tb11=373.15; %K
tb12=77.35; %K no se condensa!!!
taob1=tc1/tb1; %adimensional
taob2=tc2/tb2; %adimensional
taob3=tc3/tb3; %adimensional
taob4=tc4/tb4; %adimensional
taob5=tc5/tb5; %adimensional
taob6=tc6/tb6; %adimensional
taob7=tc7/tb7; %adimensional
taob8=tc8/tb8; %adimensional
taob9=tc9/tb9; %adimensional
taob10=tc10/tb10; %adimensional
taob11=tc11/tb11; %adimensional
g1=-5.53357241;
g2=11.0210515;
g3=-0.51243147;
g4=-10.6722729;
g5=29.4364927;
g6=-0.44101891;
m1=60.05*(1/1000);%kg/mol
m2=72.11*(1/1000);%kg/mol
m3=74.08*(1/1000);%kg/mol
m4=180.08*(1/1000);%kg/mol
m5=460.25*(1/1000);%kg/mol
m6=460.25*(1/1000);%kg/mol
m7=124.14*(1/1000);%kg/mol
m8=88.15*(1/1000);%kg/mol
m9=94.11*(1/1000);%kg/mol
m10=58.08*(1/1000);%kg/mol
m11=18.02*(1/1000);%kg/mol
%%%%EOS%%%%%
f0taob1=g1*(taob1-exp(1-taob1))+g2*(((taob1)^g3)-exp(1-taob1));
f0taob2=g1*(taob2-exp(1-taob2))+g2*(((taob2)^g3)-exp(1-taob2));
f0taob3=g1*(taob3-exp(1-taob3))+g2*(((taob3)^g3)-exp(1-taob3));
f0taob4=g1*(taob4-exp(1-taob4))+g2*(((taob4)^g3)-exp(1-taob4));
f0taob5=g1*(taob5-exp(1-taob5))+g2*(((taob5)^g3)-exp(1-taob5));
f0taob6=g1*(taob6-exp(1-taob1))+g2*(((taob6)^g3)-exp(1-taob6));
f0taob7=g1*(taob7-exp(1-taob7))+g2*(((taob7)^g3)-exp(1-taob7));
f0taob8=g1*(taob8-exp(1-taob8))+g2*(((taob8)^g3)-exp(1-taob8));
f0taob9=g1*(taob9-exp(1-taob9))+g2*(((taob9)^g3)-exp(1-taob9));
f0taob10=g1*(taob10-exp(1-taob10))+g2*(((taob10)^g3)-exp(1-taob10));
f0taob11=g1*(taob11-exp(1-taob11))+g2*(((taob11)^g3)-exp(1-taob11));
f1taob1=g4*(taob1-exp(1-taob1))+g5*(((taob1)^g6)-exp(1-taob1));
f1taob2=g4*(taob2-exp(1-taob2))+g5*(((taob2)^g6)-exp(1-taob2));
f1taob3=g4*(taob3-exp(1-taob3))+g5*(((taob3)^g6)-exp(1-taob3));
f1taob4=g4*(taob4-exp(1-taob4))+g5*(((taob4)^g6)-exp(1-taob4));
f1taob5=g4*(taob5-exp(1-taob5))+g5*(((taob5)^g6)-exp(1-taob5));
f1taob6=g4*(taob6-exp(1-taob6))+g5*(((taob6)^g6)-exp(1-taob6));
f1taob7=g4*(taob7-exp(1-taob7))+g5*(((taob7)^g6)-exp(1-taob7));
f1taob8=g4*(taob8-exp(1-taob8))+g5*(((taob8)^g6)-exp(1-taob8));
f1taob9=g4*(taob9-exp(1-taob9))+g5*(((taob9)^g6)-exp(1-taob9));
f1taob10=g4*(taob10-exp(1-taob10))+g5*(((taob10)^g6)-exp(1-taob10));
f1taob11=g4*(taob11-exp(1-taob11))+g5*(((taob11)^g6)-exp(1-taob11));
%%%estimaci�n del factor acentrico%%%
w1=(0.013162987-log(pc1)-f0taob1)/f1taob1;
w2=(0.013162987-log(pc2)-f0taob2)/f1taob2;
w3=(0.013162987-log(pc3)-f0taob3)/f1taob3;
w4=(0.013162987-log(pc4)-f0taob4)/f1taob4;
w5=(0.013162987-log(pc5)-f0taob5)/f1taob5;
w6=(0.013162987-log(pc6)-f0taob6)/f1taob6;
w7=(0.013162987-log(pc7)-f0taob7)/f1taob7;
w8=(0.013162987-log(pc8)-f0taob8)/f1taob8;
w9=(0.013162987-log(pc9)-f0taob9)/f1taob9;
w10=(0.013162987-log(pc10)-f0taob10)/f1taob10;
w11=(0.013162987-log(pc11)-f0taob11)/f1taob11;
%%%correlaci�n de pitzer%%%
T=313.15;%K
tr1=T/tc1;%K
tr2=T/tc2;%K
tr3=T/tc3;%K
tr4=T/tc4;%K
tr5=T/tc5;%K
tr6=T/tc6;%K
tr7=T/tc7;%K
tr8=T/tc8;%K
tr9=T/tc9;%K
tr10=T/tc10;%K
tr11=T/tc11;%K
%%%calores latentes de vaporizaci�n%%%
hv1=(7.08*((1-tr1)^0.354)+10.95*w1*(1-tr1)^0.456)*R*tc1*m1;%J/kg  
hv2=(7.08*((1-tr2)^0.354)+10.95*w2*(1-tr2)^0.456)*R*tc2*m2;%J/kg   
hv3=(7.08*((1-tr3)^0.354)+10.95*w3*(1-tr3)^0.456)*R*tc3*m3;%J/kg  
hv4=(7.08*((1-tr4)^0.354)+10.95*w4*(1-tr4)^0.456)*R*tc4*m4;%J/kg  
hv5=(7.08*((1-tr5)^0.354)+10.95*w5*(1-tr5)^0.456)*R*tc5*m5;%J/kg  
hv6=(7.08*((1-tr6)^0.354)+10.95*w6*(1-tr6)^0.456)*R*tc6*m6;%J/kg 
hv7=(7.08*((1-tr7)^0.354)+10.95*w7*(1-tr7)^0.456)*R*tc7*m7;%J/kg  
hv8=(7.08*((1-tr8)^0.354)+10.95*w8*(1-tr8)^0.456)*R*tc8*m8;%J/kg  
hv9=(7.08*((1-tr9)^0.354)+10.95*w9*(1-tr9)^0.456)*R*tc9*m9;%J/kg  
hv10=(7.08*((1-tr10)^0.354)+10.95*w10*(1-tr10)^0.456)*R*tc10*m10;%J/kg  
hv11=(7.08*((1-tr11)^0.354)+10.95*w11*(1-tr11)^0.456)*R*tc11*m11;%J/kg  
a1=195.75;
a2=245.97;
a3=229.04;
a4=527.97;
a5=326.70;
a6=531.25;
a7=202.39;
a8=-158.75;
a9=240.37;
a10=164.92;
a11=1779.02;
a1n=0.03298677E+02;
a2n=0.14082404E-02;
a3n=-0.03963222E-04;
a4n=0.05641515E-07;
a5n=-0.02444854E-10;
b1=3.52;
b2=4.46;
b3=3.98;
b4=3.11;
b5=2.52;
b6=3.07;
b7=4.75;
b8=4.96;
b9=4.23;
b10=4.01;
b11=0.17;
c1=-0.0015;
c2=-0.0017;
c3=-0.0015;
c4=-0.0007;
c5=-0.0010;
c6=-0.00074;
c7=-0.0019;
c8=-0.0024;
c9=-0.0017;
c10=-0.0017;
c11=0.00036;
%%%cp=a+bt+ct^2%%%
%%%% demanda energ�tica condensaci�n %%%
T1=623;%K
T2=313;%K
%%%entalp�as espec�ficas%%%
cpv1=(a1*(tb1-T1))+(b1*(1/2)*(tb1^2-T1^2))+(c1*(1/3)*(tb1^3-T1^3));%J/kg
cpl1=(a1*(T2-tb1))+(b1*(1/2)*(T2^2-tb1^2))+(c1*(1/3)*(T2^3-tb1^3));%J/kg
cpv2=(a1*(tb2-T1))+(b1*(1/2)*(tb2^2-T1^2))+(c1*(1/3)*(tb2^3-T1^3));%J/kg
cpl2=(a1*(T2-tb2))+(b1*(1/2)*(T2^2-tb2^2))+(c1*(1/3)*(T2^3-tb2^3));%J/kg
cpv3=(a1*(tb3-T1))+(b1*(1/2)*(tb3^2-T1^2))+(c1*(1/3)*(tb3^3-T1^3));%J/kg
cpl3=(a1*(T2-tb3))+(b1*(1/2)*(T2^2-tb3^2))+(c1*(1/3)*(T2^3-tb3^3));%J/kg
cpv4=(a1*(tb4-T1))+(b1*(1/2)*(tb4^2-T1^2))+(c1*(1/3)*(tb4^3-T1^3));%J/kg
cpl4=(a1*(T2-tb4))+(b1*(1/2)*(T2^2-tb4^2))+(c1*(1/3)*(T2^3-tb4^3));%J/kg
cpv5=(a1*(tb5-T1))+(b1*(1/2)*(tb5^2-T1^2))+(c1*(1/3)*(tb5^3-T1^3));%J/kg
cpl5=(a1*(T2-tb5))+(b1*(1/2)*(T2^2-tb5^2))+(c1*(1/3)*(T2^3-tb5^3));%J/kg
cpv6=(a1*(tb6-T1))+(b1*(1/2)*(tb6^2-T1^2))+(c1*(1/3)*(tb6^3-T1^3));%J/kg
cpl6=(a1*(T2-tb6))+(b1*(1/2)*(T2^2-tb6^2))+(c1*(1/3)*(T2^3-tb6^3));%J/kg
cpv7=(a1*(tb7-T1))+(b1*(1/2)*(tb7^2-T1^2))+(c1*(1/3)*(tb7^3-T1^3));%J/kg
cpl7=(a1*(T2-tb7))+(b1*(1/2)*(T2^2-tb7^2))+(c1*(1/3)*(T2^3-tb7^3));%J/kg
cpv8=(a1*(tb8-T1))+(b1*(1/2)*(tb8^2-T1^2))+(c1*(1/3)*(tb8^3-T1^3));%J/kg
cpl8=(a1*(T2-tb8))+(b1*(1/2)*(T2^2-tb8^2))+(c1*(1/3)*(T2^3-tb8^3));%J/kg
cpv9=(a1*(tb9-T1))+(b1*(1/2)*(tb9^2-T1^2))+(c1*(1/3)*(tb9^3-T1^3));%J/kg
cpl9=(a1*(T2-tb9))+(b1*(1/2)*(T2^2-tb9^2))+(c1*(1/3)*(T2^3-tb9^3));%J/kg
cpv10=(a1*(tb10-T1))+(b1*(1/2)*(tb10^2-T1^2))+(c1*(1/3)*(tb10^3-T1^3));%J/kg
cpl10=(a1*(T2-tb10))+(b1*(1/2)*(T2^2-tb10^2))+(c1*(1/3)*(T2^3-tb10^3));%J/kg
cpv11=(a1*(tb11-T1))+(b1*(1/2)*(tb11^2-T1^2))+(c1*(1/3)*(tb11^3-T1^3));%J/kg
cpl11=(a1*(T2-tb11))+(b1*(1/2)*(T2^2-tb11^2))+(c1*(1/3)*(T2^3-tb11^3));%J/kg
h1=(cpv1-hv1+cpl1)*(x1*F*p1)/1000;%kW
h2=(cpv2-hv2+cpl2)*(x2*F*p2)/1000;%kW
h3=(cpv3-hv3+cpl3)*(x3*F*p3)/1000;%kW
h4=(cpv4-hv4+cpl4)*(x4*F*p4)/1000;%kW
h5=(cpv5-hv5+cpl5)*(x5*F*p5)/1000;%kW
h6=(cpv6-hv6+cpl6)*(x6*F*p6)/1000;%kW
h7=(cpv7-hv7+cpl7)*(x7*F*p7)/1000;%kW
h8=(cpv8-hv8+cpl8)*(x8*F*p8)/1000;%kW
h9=(cpv9-hv9+cpl9)*(x9*F*p9)/1000;%kW
h10=(cpv10-hv10+cpl10)*(x10*F*p10)/1000;%kW
h11=(cpv11-hv11+cpl11)*(x11*F*p11)/1000;%kW
h12=((a1n*(T2-T1))+(a2n*(1/2)*(T2^2-T1^2))+(a3n*(1/3)*(T2^3-T1^3))+(a4n*(1/4)*(T2^4-T1^4))+(a5n*(1/5)*(T2^5-T1^5)))*8.314*28.01*(x12*F*p12);
%%%demanda energ�tica%%%
qneto=h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12;%kW
%%condensaci�n total de los productos condensables%%
%%%%%%fluido de enfriamiento%%%%%%
%%%primera ley de la termodin�mica%%%
%%%m*cp*dt=q%%%
%%%insertar datos de la naturaleza del fluido%%%
mw=18.02*(1/1000);%kg/mol
tbw=373.15/1000;%K
thin=298.15/1000;%K
thout=473.15/1000;%K
awl=-203.60;         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bwl=1523.29;         %                                               %
cwl=-3196.41;        % Las constantes para calcular el polinomio del %
dwl=2474.45;         % cp fueron obtenidas de la base de datos de    %
ewl=3.85;            % NIST en el siguiente enlace:                  %
                     % https://webbook.nist.gov/cgi/cbook.cgi?       %
                     % ID=C7732185&Type=JANAFL&Plot=on               %                               
                     %                                               %
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hvw=1939.7*(1/1000);%J/mol
cpest=awl+(bwl*thin)+(cwl*(thin^2))+(dwl*(thin^3))+(ewl/(thin^2));
cplw=(awl*(tbw-thin)+((1/2)*(bwl*tbw^2-thin^2))+((1/3)*cwl*(tbw^3-thin^3))+(((1/4)*dwl*(tbw^4-thin^4))-(ewl*(tbw^-1-thin^-1))));%J/mol 
cpgw=(awl*(thout-tbw)+((1/2)*bwl*(thout^2-tbw^2))+((1/3)*cwl*(thout^3-tbw^3))+(((1/4)*dwl*(thout^4-tbw^4))-(ewl*(thout^-1-tbw^-1))));%J/mol
hw=((cplw+hvw+cpgw)*(1/1000))/mw;%kJ/kg
mfw=-qneto/hw;%kg/s
volw=mfw/992.3; %m3/s