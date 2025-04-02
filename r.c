#include <cstdio>//g++ r.c -o r-lm;./r
#include <math.h>
#include <complex>
#include <vector>
#include <algorithm>
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
V N(V v){F l=L(v);K R(v,w(1./l));}
struct E{I f;F a;F p;B<F>c;};
F h=.5;
F y3=.01;
F e2=.001;
V U(F t){F x=cos(2.*P*(t));F y=cos(2.*P*(t+.33));F z=cos(2.*P*(t+.67));K A(w(h),R(w(h),w(x,y,z)));}
J<B<F>>q={{-h,h},{-h,-h},{0.,-.85},{h,-h},{h,h},{-h,h},{h,-h},{-h,-h},{h,h},{-h,h}};
J<B<F>>T(J<B<F>>&p,I b){F t=0;J<F>s;O(I i=0;i<9;++i){F l=abs(p[i+1]-p[i]);s.PB(l);t+=l;}J<B<F>>k;F step=t/b,d=0;I g=0;O(I i=0;i<b;++i,d+=step){while(g<8&d>s[g])d-=s[g++];F f=d/s[g];k.PB(p[g]+(p[g+1]-p[g])*f);}K k;}
J<B<F>>X=T(q,64);J<E>Q;
void g(){I N=64;Q.clear();O(I k=0;k<N;++k){B<F>s(0.,0.);O(I n=0;n<N;++n){F h=2*P*k*n/N;s+=X[n]*exp(B<F>(0,-h));}s/=(F)N;I f=(k<=N/2)?k:k-N;Q.PB({f,abs(s),atan2(s.imag(),s.real()),s});}sort(Q.begin(),Q.end(),[](E a,E b){K a.a>b.a;});}
F tl[W][W]={0};
F Z(V p,V c,F r){K L(A(p,{-c.x,-c.y,-c.z}))-r;}
F Y(V p,J<V>&c,V t,F r,I *o,I *j){F m=1e9;I b=-1;I id=b;O(I i=0;i<64;++i){F f=max(Q[i].a,y3);F d=Z(p,c[i],f);if(d<m){m=d;id=1;b=i;}}F d=Z(p,t,r);if(d<m){m=d;id=2;b=-1;}*o=id;*j=b;K m;}
V l(F t,J<V>&c,F &x){B<F>p(0.,0.);c.clear();F r=0.;
O(auto &e:Q)r+=e.a;F z=r+.1,d=r/64;O(I i=0;i<64;++i){E &e=Q[i];F a=e.f*t+e.p;B<F>o=e.a*exp(B<F>(0,a));c.PB(w(p.real(),p.imag(),z));p+=o;z-=d;if(i==63)x=max(e.a,y3);}K w(p.real(),p.imag(),z);}
V C(V p,J<V>&c,V t,F r){I d,b;auto m=[&](F x,F y,F z){K Y(A(p,w(x,y,z)),c,t,r,&d,&b);};K N(w(m(e2,0,0)-m(-e2,0,0),m(0,e2,0)-m(0,-e2,0),m(0,0,e2)-m(0,0,-e2)));}
main(){
g();I f,z,s,a,c=400,x,y,i,e,h;F t,u,v,d,k,tr,bt,at,r,g,b,D,br;
J<V>CT;V ti,ld,n,m,p,j=w(0,0,-1),rd;
char G[16];FILE*o;
O(f=0;f<c;f++){
sprintf(G,"%03d.ppm",f);
o=fopen(G,"wb");
fputs("P6\n512 512\n255\n",o);
t=2*P*f/c;ti=l(t,CT,tr);D=3;t*=1.5;
ld=N(w(cos(t)*D,sin(t*h)*1.5,sin(t)*D));
O(y=0;y<W;y++)O(x=0;x<W;x++){
u=(2.*x-W)/W;v=(2.*y-W)/W;rd=N(w(u,v,1));d=0;a=0;
O(i=c;d<3&&i--;d+=k){
p=A(j,R(rd,w(d)));
if((k=Y(p,CT,ti,tr,&s,&z))<e2){a=s;break;}
}
r=g=b=0;
if(a){
p=A(j,R(rd,w(d)));n=C(p,CT,ti,tr);
br=.3+fmax(0.,n.x*ld.x+n.y*ld.y+n.z*ld.z)*.7;
bt=(z>=0)?(F)z/(63):0;
at=fmod((a==1?sqrt(bt):1)+f*y3,1);m=U(at);
r=br*m.x;g=br*m.y;b=br*m.z;
if(a-1){e=(u+1)*256;h=(v+1)*256;if(e>=0&e<W&h>=0&h<W)tl[h][e]=1;}}
k=tl[y][x];
r+=k;r=r>1?1:r;g*=1-k;g=g>1?1:g;b*=1-k;b=b>1?1:b;
fputc(255*r,o);fputc(255*g,o);fputc(255*b,o);}}}