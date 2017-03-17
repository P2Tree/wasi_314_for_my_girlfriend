function E=fun(rs, ri, Ci)
    q = 500000;
    k = 1000;
    xi = ri(1, :);
    yi = ri(2, :);
    V = ri(3, :);
    theta = ri(4, :);
    distance = sqrt((rs(1)-xi).^2+(rs(2)-yi).^2);        
    deltax = (rs(1) - xi).*cos(theta) + (rs(2) - yi).*sin(theta);     
    ci=q/(2*pi*k)*1./(distance).*exp(-V./(2*k).*(distance)-deltax);
    E = sum((Ci - ci).^2);
        