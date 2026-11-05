import numpy as np; from PIL import Image, ImageDraw
im=np.array(Image.open('19930086559/pages/p14-14.png').convert('L')).astype(int)
grid=np.array([952,1014,1076,1138,1199,1261,1322,1384,1445,1507,1569,1611,1631,1693,1754,1815,1878,1939,2001,2062,2124,2186,2248,2310,2372,2434])
ax=[920,1000,1100,1200,1290,1400,1500,1600,1700,1800]; ay=[1109,1252,1415,1569,1690,1834,1946,2047,2135,2212]
pts=[]
for x in range(906,1851,4):
    yp=np.interp(x,ax,ay)
    col=im[int(yp-30):int(yp+31),x]; ys=np.where(col<130)[0]+int(yp-30)
    g=[]
    for y in ys:
        if g and y-g[-1][-1]<=2: g[-1].append(y)
        else: g.append([y])
    c=[np.mean(a) for a in g if len(a)<=8 and np.min(abs(grid-np.mean(a)))>4]
    if c:
        y=min(c,key=lambda v:abs(v-yp))
        if abs(y-yp)<20: pts.append((x,y))
pts=np.array(pts)
# smooth outliers: residual vs local quadratic
keep=[]
for i,(x,y) in enumerate(pts):
    m=(abs(pts[:,0]-x)<60)&(abs(pts[:,0]-x)>0)
    p=np.polyfit(pts[m,0],pts[m,1],2); keep.append(abs(np.polyval(p,x)-y)<4)
pts=pts[np.array(keep)]
np.save('beta_pts.npy',pts)
o=Image.open('19930086559/pages/p14-14.png').convert('RGB'); d=ImageDraw.Draw(o)
for x,y in pts: d.ellipse([x-3,y-3,x+3,y+3],outline=(255,0,0))
o.crop((880,900,1880,2460)).resize((700,1092)).save('ov51_beta.png')
print(len(pts))
