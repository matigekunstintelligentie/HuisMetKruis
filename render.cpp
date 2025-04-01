#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
// g++ -O2 render.cpp -o render -lm;./render
using namespace std;
#define W 512
#define FRAMES 400
#define F float
#define I int
#define B complex
#define J vector
#define K return
#define O for
F P=3.141592653589793f;
typedef struct{F x,y,z;} V;
V w(F x,F y,F z){K{x,y,z};}
V w(F x){K{x,x,x};}
V A(V a,V b){K w(a.x+b.x,a.y+b.y,a.z+b.z);}
V M(V a,F s){K w(a.x*s,a.y*s,a.z*s);}
V R(V a,V b){K w(a.x*b.x,a.y*b.y,a.z*b.z);}
F D(V a,V b){K a.x*b.x+a.y*b.y+a.z*b.z;}
F L(V v){K sqrtf(D(v,v));}
V N(V v){F l=L(v);K(l>0)?M(v,1.0f/l):v;}
struct E{I f;F a;F p;B<F>c;};
F h=.5;
V U(F t){F x=cosf(2.*P*(t));F y=cosf(2.*P*(t+.33));F z=cosf(2.*P*(t+.67));K A(w(h),R(w(h),w(x,y,z)));}
J<B<F>>q={{-h,h},{-h,-h},{0.,-.85},{h,-h},{h,h},{-h,h},{h,-h},{-h,-h},{h,h},{-h,h}};
J<B<F>>T(J<B<F>>& p,I b){F t=0.;J<F>s;O(I i=0;i<p.size()-1;++i){F l=abs(p[i+1]-p[i]);s.push_back(l);t+=l;}J<B<F>>k;F step=t/b;F c=0.;I g=0;F q=c;O(I i=0;i<b;++i){while(g<s.size()&&q+s[g]<c){q+=s[g];g++;}if(g>=s.size())break;F l=(c-q)/s[g];B<F>a=p[g];B<F>b=p[g+1];B<F>x=a+(b-a)*l;k.push_back(x);c+=step;}K k;}
J<B<F>>X;
J<E>Q;
void g(){I N=X.size();Q.clear();O(I k=0;k<N;++k){B<F>s(0.,0.);O(I n=0;n<N;++n){F h=2*P*k*n/N;s+=X[n]*exp(B<F>(0,-h));}s/=static_cast<F>(N);I f=(k<=N/2)?k:k-N;Q.push_back({f,abs(s),atan2(s.imag(),s.real()),s});}sort(Q.begin(),Q.end(),[](E a,E b){K a.a>b.a;});}
F trail[W][W]={0};
F Z(V p,V c,F r){K L(A(p,{-c.x,-c.y,-c.z}))-r;}
F Y(V p,J<V>&c,V t,F r,I *o,I *j){F m=1e9;I id=-1;I b=-1;O(I i=0;i<c.size();++i){F f=fmaxf(Q[i].a,.01);F d=Z(p,c[i],f);if(d<m){m=d;id=1;b=i;}}F d=Z(p,t,r);if(d<m){m=d;id=2;b=-1;}*o=id;*j=b;K m;}
V l(F t,J<V>&c,F &x){B<F>p(0.,0.);c.clear();F r=0.;
O(auto &e:Q)r+=e.a;
F z=r+.1;
F d=r/Q.size();
O(I i=0;
i<Q.size();
++i){
E &e=Q[i];
F a=e.f*t+e.p;
B<F>o=e.a*exp(B<F>(0,a));
c.push_back(w(p.real(),p.imag(),z));
p+=o;
z-=d;
if(i==Q.size()-1)x=fmaxf(e.a,.01);
}
K w(p.real(),p.imag(),z);
}
V C(V p,J<V>&c,V t,F r){F e=.001;I d,b;F x=Y(A(p,w(e,0,0)),c,t,r,&d,&b)-Y(A(p,w(-e,0,0)),c,t,r,&d,&b);F y=Y(A(p,w(0,e,0)),c,t,r,&d,&b)-Y(A(p,w(0,-e,0)),c,t,r,&d,&b);F z=Y(A(p,w(0,0,e)),c,t,r,&d,&b)-Y(A(p,w(0,0,-e)),c,t,r,&d,&b);K N(w(x,y,z));}

void write_ppm_header(FILE *f){
    fprintf(f,"P6\n%d %d\n255\n",W,W);
}

V get_light(F t){
    F angle=t*1.5;
    F radius=3.;
    F x=cosf(angle)*radius;
    F y=sinf(angle*h)*1.5;
    F z=sinf(angle)*radius;
    K N(w(x,y,z));
}

I main(){
    X=T(q,64);
    g();
    O(I frame=0;frame<FRAMES;++frame){
        char fname[64];
        snprintf(fname,sizeof(fname),"frame_%03d.ppm",frame);
        FILE *out=fopen(fname,"wb");
        if(!out) continue;
        write_ppm_header(out);

        F t=2*P*frame/FRAMES;
        F tip_radius;
        J<V>centers;
        V tip=l(t,centers,tip_radius);

        O(I y=0;y<W;y++)
            O(I x=0;x<W;x++)
                trail[y][x] *=.999;

        V light_dir=get_light(t);

        O(I y=0;y<W;y++){
            O(I x=0;x<W;x++){
                I shape_id=0;
                I shape_index=-1;
    
                F u=(2.*x-W)/W;
                F v=(2.*y-W)/W;
                V ro=w(0,0,-1);
                V rd=N(w(u,v,1));
                F d=0.;
                I id=0;

                O(I i=0;i<100;i++){
                    V p=A(ro,M(rd,d));
                    F dist=Y(p,centers,tip,tip_radius,&shape_id,&shape_index);

                    if(dist<.001){ id=shape_id;break;}
                    if(d>3.) break;
                    d+=dist;
                }

                unsigned char r=0,g=0,b=0;
                if(id==1 || id==2){
                    V p=A(ro,M(rd,d));
                    V n=C(p,centers,tip,tip_radius);
                    F ambient=.3;
                    F brightness=ambient+fmaxf(0.,D(n,light_dir))*.7;

                    F ball_t=(shape_index>=0) ? (F)shape_index/(F)(centers.size()-1):0.;

                    if(id==1){
                        F animated_t=fmodf(sqrt(ball_t)+frame*.01,1.);
                        V color=U(animated_t);
                        r=(unsigned char)(brightness*255*color.x);
                        g=(unsigned char)(brightness*255*color.y);
                        b=(unsigned char)(brightness*255*color.z);
                    } else if(id==2){
                        F tip_t=fmodf(sqrt(1.)+frame*.01,1.); // which is just fmodf(1.0+frame*0.01f,1.0f)

                        V color=U(tip_t);
                        r=(unsigned char)(brightness*255*color.x);
                        g=(unsigned char)(brightness*255*color.y);
                        b=(unsigned char)(brightness*255*color.z);
                        I tx=(int)((u+1)*h*W);
                        I ty=(int)((v+1)*h*W);
                        if(tx>=0 && tx<W && ty>=0 && ty<W)
                            trail[ty][tx]=1.;
                    }
                }

                F trail_strength=trail[y][x];
                if(trail_strength>.001){
                    // Mix original color with red trail
                    F tr=trail_strength*255.;
                    r=(unsigned char)clamp((F)r+tr,0.f,255.f);
                    g=(unsigned char)clamp((F)g*(1.-trail_strength),0.,255.);
                    b=(unsigned char)clamp((F)b*(1.-trail_strength),0.,255.);
                }

                fputc(r,out);fputc(g,out);fputc(b,out);
            }
        }

        fclose(out);
    }

    K 0;
}