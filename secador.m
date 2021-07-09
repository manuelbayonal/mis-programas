%secador de biomasa%
%flujo de biomasa seca%
f2=4.7222;%kg/s
%porcentaje de humedad
wcsi=0.25;%inicial
wcsf=0.08; %final
%f1+f2=f3+f4
%temperatura de entrada del aire
Tairin=150;%�C
Tairout=120;%�C
wcai=0.55/(1+0.55);%inicial %carta psicrom�trica en linea: http://www.flycarpet.net/en/PsyOnline.
wcaf=0.9/(1+0.9);%final %carta psicrom�trica en l�nea: http://www.flycarpet.net/en/PsyOnline.
%f1+f3=f2+f4
%balances por componente
%biomasa
f1=f2*(1-wcsf)/(1-wcsi);%kg/s
%f1+f3-f2=f4
%agua
f3=(wcsi*f1-wcsf*f2+wcaf*(f2-f1))/(wcaf-wcai);%kg/s
f4=f1+f3-f2;%kg/s
