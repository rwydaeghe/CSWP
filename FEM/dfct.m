function d = dfct(p,f)
%y=x^2--> d=y-f(x) is positief als boven de functie, dus negatief als onder. Neem dus dit
%A min B: max(dA,-dB). En B=negatieve y halfvlak--> dB=y
%(A min B) min C, met C = negatieve x halfvlak
y=p(:,2);
x=p(:,1);
d=max(max(y-f(x),-y),-x);


