#include <cstdio>
#include <math.h>
#include <complex>
#include <vector>
#include <algorithm>
//g++ r.c -o r-lm;./r
using namespace std;
#define F float
#define I int
#define B complex
#define J vector
#define K return
#define O for
#define PB push_back
const I W=512;
F P=3.1416;
typedef struct{F x,y,z;} V;
V w(F x,F y,F z){K{x,y,z};}
V w(F x){K{x,x,x};}
V A(V a,V b){K w(a.x+b.x,a.y+b.y,a.z+b.z);}
V R(V a,V b){K w(a.x*b.x,a.y*b.y,a.z*b.z);}
F L(V v){K sqrt(v.x*v.x+v.y*v.y+v.z*v.z);}
V N(V v){F l=L(v);K(l>0)?R(v,w(1.0f/l)):v;}
struct E{I f;F a;F p;B<F>c;};
F h=.5;
F y3=.01;
V U(F t){F x=cos(2.*P*(t));F y=cos(2.*P*(t+.33));F z=cos(2.*P*(t+.67));K A(w(h),R(w(h),w(x,y,z)));}
J<B<F>>q={{-h,h},{-h,-h},{0.,-.85},{h,-h},{h,h},{-h,h},{h,-h},{-h,-h},{h,h},{-h,h}};
J<B<F>>T(J<B<F>>&p,I b){F t=0;J<F>s;O(I i=0;i+1<p.size();++i){F l=abs(p[i+1]-p[i]);s.PB(l);t+=l;}J<B<F>>k;F step=t/b,d=0;I g=0;O(I i=0;i<b;++i,d+=step){while(g+1<s.size()&d>s[g])d-=s[g++];F f=d/s[g];k.PB(p[g]+(p[g+1]-p[g])*f);}K k;}
J<B<F>>X;J<E>Q;
void g(){I N=64;Q.clear();O(I k=0;k<N;++k){B<F>s(0.,0.);O(I n=0;n<N;++n){F h=2*P*k*n/N;s+=X[n]*exp(B<F>(0,-h));}s/=(F)N;I f=(k<=N/2)?k:k-N;Q.PB({f,abs(s),atan2(s.imag(),s.real()),s});}sort(Q.begin(),Q.end(),[](E a,E b){K a.a>b.a;});}
F tl[W][W]={0};
F Z(V p,V c,F r){K L(A(p,{-c.x,-c.y,-c.z}))-r;}
F Y(V p,J<V>&c,V t,F r,I *o,I *j){F m=1e9;I id=-1;I b=-1;O(I i=0;i<64;++i){F f=max(Q[i].a,y3);F d=Z(p,c[i],f);if(d<m){m=d;id=1;b=i;}}F d=Z(p,t,r);if(d<m){m=d;id=2;b=-1;}*o=id;*j=b;K m;}
V l(F t,J<V>&c,F &x){B<F>p(0.,0.);c.clear();F r=0.;
O(auto &e:Q)r+=e.a;F z=r+.1;F d=r/64;O(I i=0;i<64;++i){E &e=Q[i];F a=e.f*t+e.p;B<F>o=e.a*exp(B<F>(0,a));c.PB(w(p.real(),p.imag(),z));p+=o;z-=d;if(i==Q.size()-1)x=max(e.a,y3);}K w(p.real(),p.imag(),z);}
V C(V p,J<V>&c,V t,F r){F e=.001;I d,b;auto m=[&](F x,F y,F z){K Y(A(p,w(x,y,z)),c,t,r,&d,&b);};K N(w(m(e,0,0)-m(-e,0,0),m(0,e,0)-m(0,-e,0),m(0,0,e)-m(0,0,-e)));}
I main(){
X=T(q,64);g();I f,z,s,id,fr=400,x,y,i,tx,ty;F t,u,v,d,k,tr,bt,at,r,g,b,rad,br;
J<V>CT;V ti,ld,n,m,p,ro=w(0,0,-1),rd;
char fn[16];FILE*o;
O(f=0;f<fr;f++){
snprintf(fn,16,"%03d.ppm",f);
o=fopen(fn,"wb");
fprintf(o,"P6\n%d %d\n255\n",W,W);
t=2*P*f/fr;ti=l(t,CT,tr);rad=3;t*=1.5;
ld=N(w(cos(t)*rad,sin(t*h)*1.5,sin(t)*rad));
O(y=0;y<W;y++)O(x=0;x<W;x++){
u=(2.*x-W)/W;v=(2.*y-W)/W;rd=N(w(u,v,1));d=0;id=0;
O(i=fr;i--;){
p=A(ro,R(rd,w(d)));
if((k=Y(p,CT,ti,tr,&s,&z))<.001){id=s;break;}
if(d>3)break;d+=k;}
r=g=b=0;
if(id){
p=A(ro,R(rd,w(d)));n=C(p,CT,ti,tr);
br=.3+fmax(0.f,n.x*ld.x+n.y*ld.y+n.z*ld.z)*.7;
bt=(z>=0)?(F)z/(CT.size()-1):0;
at=fmod((id==1?sqrt(bt):1)+f*y3,1);m=U(at);
r=br*m.x;g=br*m.y;b=br*m.z;
if(id-1){tx=(u+1)*256;ty=(v+1)*256;if(tx>=0&tx<W&ty>=0&ty<W)tl[ty][tx]=1;}}
k=tl[y][x];
r=fminf(r+k,1.f);g=fminf(g*(1-k),1.f);b=fminf(b*(1-k),1.f);
fputc((I)(255*r),o);fputc((I)(255*g),o);fputc((I)(255*b),o);}}fclose(o);}