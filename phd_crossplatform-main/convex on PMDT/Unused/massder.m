function mdot = massder(t,m,us, param)

ux = interp1(param.tvec, us(:,1), t);
uy = interp1(param.tvec, us(:,2), t);
uz = interp1(param.tvec, us(:,3), t);

u = [ux;uy;uz];
unorm = norm(u);

mdot = -m*unorm/param.Isp/param.g0;
