function [V,F,b,bc] = airfoil(n)
%AIRFOIL Construct an arfoil with resolution n
%fl is set to false if this is a normal domain, true if it is a domain with
%a single hole in the middle.

if nargin<1
    n = 100;
end

hlead=0.05; htrail=0.04; hmax=2*pi/n; circx=2; circr=4;
a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];

fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);

fixx=1-htrail*cumsum(1.3.^(0:4)');
fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
box=[circx-circr,-circr; circx+circr,circr];
h0=min([hlead,htrail,hmax]);

[p,t]=distmesh2d(fd,fh,h0,box,fix);

wb = unique(outline(t));

outerb = normrow(p(wb,:) - [2 0]) > 3;
innerb_lo = ~outerb & (p(wb,2)<0);
innerb_hi = ~outerb & (p(wb,2)>=0);

id_outerb = wb(outerb);
id_innerb_lo = wb(innerb_lo);
id_innerb_hi = wb(innerb_hi);

V = {p};
F = {t};
b = {[id_outerb; id_innerb_lo; id_innerb_hi]};
bc = {[repmat([255,255,191]/255,numel(id_outerb),1); ...
    repmat([215,48,39]/255,numel(id_innerb_lo),1); ...
    repmat([69,117,180]/255,numel(id_innerb_hi),1)]};

end

