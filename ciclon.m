%Dimensi�n (L*W*H):800x800x1633mm
%fuente de las dimensiones: https://spanish.alibaba.com/product-detail/dust-filter-filtrator-with-small-cyclone-separator-industrial-dust-collector-fume-extraction-system-60606308109.html?spm=a2700.galleryofferlist.normal_offer.d_image.5e754e9bRbVz3v&s=p
L=0.8;%m
Z=0.8;%m
H=1.63;%m
Ne=(1/H)*(L+(Z/2));%revoluciones efectivas.
%Punto de corte
dpc=350;%micrometros
%teniendo en cuenta que, el modelo de koufopanos no considera la
%contracci�n de las part�culas durante la pir�lisis:
dpj=85000;
%eficiencia fraccional
niuj=1/(1+(dpc/dpj)^2);
%% separaci�n total de la fase s�lida debido al punto de corte