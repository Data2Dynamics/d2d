function error_ellipse3(C, mu, color, conf)
xfc = linspace(mu(1)-10*sqrt(C(1,1)), mu(1)+10*sqrt(C(1,1)), 100);
yfc = linspace(mu(2)-10*sqrt(C(2,2)), mu(1)+10*sqrt(C(2,2)), 100);
zfc = linspace(mu(3)-10*sqrt(C(3,3)), mu(1)+10*sqrt(C(3,3)), 100);
[XFC, YFC, ZFC] = meshgrid(xfc, yfc, zfc);
P = mvnpdf([XFC(:) YFC(:) ZFC(:)], mu, C);
P = P / mvnpdf(mu, mu, C);
P = reshape(P, size(XFC));
p = patch(isosurface(XFC, YFC, ZFC, P, 1-conf));
isonormals(XFC, YFC, ZFC, P, p)
set(p,'FaceColor',color,'EdgeColor','none');
set(p,'FaceAlpha',0.3);
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud