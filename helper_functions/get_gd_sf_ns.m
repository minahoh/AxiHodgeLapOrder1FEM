function [gd,sf,ns] = get_gd_sf_ns(xs, ys)


if xs == [0,1,1,0] & ys == [0,0,1,1]
    gd = [2;4;0;1;1;0;0;0;1;1];
    sf = 'P1';
    ns = [80;49];
end
