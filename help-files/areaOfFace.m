function A = areaOfFace(x,y,z)
    x = x(:)';
    y = y(:)';
    z = z(:)';
    ons = [1 1 1];
    A = 0.5*sqrt(det([x;y;ons])^2 + det([y;z;ons])^2 + det([z;x;ons])^2);
end