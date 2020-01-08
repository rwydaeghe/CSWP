function step_SIT_SIP_impedancepml(nx,ny,c,dx,dy,dt,npml,nsf,it)
global ox oy p pp
global qx qy q
global oxref pref
Z=c; %surface impedance for non-reflecting boundary

km = 2000;
n = nx+2*npml+2*nsf;

for i=1:nx+2*npml+2*nsf
    for j=1:ny+2*npml+2*nsf
        if i <= npml
            kx = km*((npml-i+1/2)/npml)^4;
        elseif i >= nx+2*nsf+npml+1
            kx = km*((i-nx-2*nsf-npml-1/2)/npml)^4;
        else
            kx = 0;
        end
        if j <= npml
            ky = km*((npml-j+1/2)/npml)^4;
        elseif j >= ny+npml+2*nsf+1
            ky = km*((j-ny-npml-2*nsf-1/2)/npml)^4;
        else
            ky = 0;
        end
        if i==npml+nsf
            if (j > npml+nsf) && (j<=n-npml-nsf)
                pp(j,i) = ((1-kx*dt/2)*pp(j,i) - c^2*dt/dy*(ox(j,i+1)-oxref(it,i+1)-ox(j,i)))/(1+kx*dt/2);
            else
                pp(j,i) = ((1-kx*dt/2)*pp(j,i) - c^2*dt/dy*(ox(j,i+1)-ox(j,i)))/(1+kx*dt/2);
            end
        elseif i == n-npml-nsf+1
            if (j > npml+nsf) && (j<=n-npml-nsf)
                pp(j,i) = ((1-kx*dt/2)*pp(j,i) - c^2*dt/dy*(ox(j,i+1)+oxref(it,i)-ox(j,i)))/(1+kx*dt/2);
            else
                pp(j,i) = ((1-kx*dt/2)*pp(j,i) - c^2*dt/dy*(ox(j,i+1)-ox(j,i)))/(1+kx*dt/2);
            end
        else
            pp(j,i) = ((1-kx*dt/2)*pp(j,i) - c^2*dt/dy*(ox(j,i+1)-ox(j,i)))/(1+kx*dt/2);
        end
     
        %p(j,i) = ((1-ky*dt/2)*p(j,i) - c^2*dt/dy*(oy(j+1,i)-oy(j,i)))/(1+ky*dt/2);
        
    end
end

for i=2:nx+2*npml+2*nsf
    for j=1:ny+2*npml+2*nsf
        if i <= npml
            kx = km*((npml-i)/npml)^4;
        elseif i >= nx+npml+2*nsf+1
            kx = km*((i-nx-2*nsf-npml)/npml)^4;
        else
            kx = 0;
        end
        if i == npml+nsf+1
            if (j > npml+nsf) && (j<=n-npml-nsf)
                ox(j,i) = ((1-kx*dt/2)*ox(j,i) - dt*(1/dx*(p(j,i)+pp(j,i)-p(j,i-1)-pref(it+1,i-1)-pp(j,i-1))))/(1+kx*dt/2); 
            else
                ox(j,i) = ((1-kx*dt/2)*ox(j,i) - dt*(1/dx*(p(j,i)+pp(j,i)-p(j,i-1)-pp(j,i-1))))/(1+kx*dt/2); 
            end
        elseif i == n+1-npml-nsf
            if (j > npml+nsf) && (j<=n-npml-nsf)
                ox(j,i) = ((1-kx*dt/2)*ox(j,i) - dt*(1/dx*(p(j,i)+pref(it+1,i)+pp(j,i)-p(j,i-1)-pp(j,i-1))))/(1+kx*dt/2);
            else 
                ox(j,i) = ((1-kx*dt/2)*ox(j,i) - dt*(1/dx*(p(j,i)+pp(j,i)-p(j,i-1)-pp(j,i-1))))/(1+kx*dt/2); 
            end
        else
            ox(j,i) = ((1-kx*dt/2)*ox(j,i) - dt*(1/dx*(p(j,i)+pp(j,i)-p(j,i-1)-pp(j,i-1))))/(1+kx*dt/2);
        end
    end
end

for j=1:ny+2*npml+2*nsf
    ox(j,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,1)-2*dt/dx*(p(j,1)+pp(j,1)));
    ox(j,nx+2*npml+2*nsf+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,nx+2*npml+2*nsf+1)+2*dt/dx*(p(j,nx+2*npml+2*nsf)+pp(j,nx+2*npml+2*nsf)));
end


% for i=1:nx+2*npml+2*nsf
%     for j = 2:ny+2*npml+2*nsf
%         if j <= npml
%             ky = km*((npml-j)/npml)^4;
%         elseif j >= ny+npml+2*nsf+1
%             ky = km*((j-ny-2*nsf-npml)/npml)^4;
%         else
%             ky = 0;
%         end
%         oy(j,i) = ((1-ky*dt/2)*oy(j,i) - dt*(1/dy*(p(j,i)+pp(j,i)-p(j-1,i)-pp(j-1,i))))/(1+ky*dt/2);
%     end
% end
% 
% for i=1:nx+2*npml+2*nsf
%     oy(1,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(1,i)-2*dt/dx*(p(1,i)+pp(1,i)));
%     oy(ny+2*npml+2*nsf+1,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(ny+2*npml+2*nsf+1,i)+2*dt/dx*(p(ny+2*npml+2*nsf,i)+pp(ny+2*npml+2*nsf,i)));
% end


oy = oy.*qy;
ox = ox.*qx;
p = p.*q;
pp = pp.*q;


end
        








%update ox oy p

 