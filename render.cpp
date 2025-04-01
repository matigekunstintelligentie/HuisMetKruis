#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

// g++ -O2 render.cpp -o render -lm 
// ./render
using namespace std;

#define W 512
#define H 512
#define FRAMES 400
#define F float
#define I int
#define B complex
#define J vector

F P = 3.141592653589793f;

typedef struct {F x,y,z;} V;

V w(F x,F y,F z){return{x,y,z};}
V w(F x){return{x,x,x};}
V A(V a,V b){return w(a.x+b.x, a.y + b.y, a.z + b.z);}
V M(V a,F s){return w(a.x*s,a.y*s,a.z*s);}
V R(V a,V b){return w(a.x*b.x,a.y*b.y,a.z*b.z);}
F D(V a,V b){return a.x*b.x + a.y*b.y + a.z*b.z;}
F L(V v){return sqrtf(D(v,v));}
V N(V v){F l=L(v);return(l>0)?M(v,1.0f/l):v;}

struct E{I f;F a;F p;B<F> c;};

V U(F t){F x=cosf(2.*P*(t));F y=cosf(2.*P*(t+.33));F z=cosf(2.*P*(t+.67));return A(w(.5),R(w(.5),w(x,y,z)));}

F h=.5;
J<B<F>> q={{-h,h},{-h,-h},{0.,-.85},{h,-h},{h,h},{-h,h},{h,-h},{-h,-h},{h,h},{-h,h}};

J<B<F>> T(J<B<F>>& p, I b) {F t=0.;J<F>s;for(I i=0;i<p.size()-1;++i){F l=abs(p[i+1]-p[i]);s.push_back(l);t+=l;}J<B<F>>k;F step=t/b;F c=0.;I g=0;F q=c;for(I i=0;i<b;++i){while(g<s.size()&&q+s[g]<c){q+=s[g];g++;}if(g>=s.size())break;F l=(c-q)/s[g];B<F>a=p[g];B<F>b=p[g+1];B<F>x=a+(b-a)*l;k.push_back(x);c+=step;}return k;}

J<B<F>> path;
J<E> epicycles;

void compute_dft() {
    I N = path.size();
    epicycles.clear();
    for (I k = 0; k < N; ++k) {
        B<F> sum(0., 0.);
        for (I n = 0; n < N; ++n) {
            F phi = 2 * P * k * n / N;
            sum += path[n] * exp(B<F>(0, -phi));
        }
        sum /= static_cast<F>(N);
        I f = (k <= N/2) ? k : k - N;
        epicycles.push_back({f, abs(sum), atan2(sum.imag(), sum.real()), sum});
    }
    sort(epicycles.begin(), epicycles.end(), [](E a, E b) {
        return a.a > b.a;
    });
}

F trail[H][W] = {0};

F sdSphere(V p, V center, F r) {
    return L(A(p, {-center.x,-center.y,-center.z})) - r;
}

F map(V p, J<V> &centers, V tip, F tip_radius, I *out_id, I *out_index) {
    F min_dist = 1e9;
    I id = -1;
    I closest_index = -1;

    for (I i = 0; i < centers.size(); ++i) {
        F radius = fmaxf(epicycles[i].a, .01);
        F d = sdSphere(p, centers[i], radius);
        if (d < min_dist) {
            min_dist = d;
            id = 1;
            closest_index = i;
        }
    }
    F d_tip = sdSphere(p, tip, tip_radius);
    if (d_tip < min_dist) {
        min_dist = d_tip;
        id = 2;
        closest_index = -1;
    }
    *out_id = id;
    *out_index = closest_index;
    return min_dist;
}

V get_tip_and_centers(F t, J<V> &centers_out, F &tip_radius_out) {
    B<F> pos(0., 0.);
    centers_out.clear();

    F total_radius = 0.;
    for (auto &e : epicycles)
        total_radius += e.a;

    F z = total_radius + .1;
    F dz = total_radius / epicycles.size();

    for (I i = 0; i < epicycles.size(); ++i) {
        E &e = epicycles[i];
        F angle = e.f * t + e.p;
        B<F> offset = e.a * exp(B<F>(0, angle));
        centers_out.push_back(w(pos.real(), pos.imag(), z));
        pos += offset;
        z -= dz;
        if (i == epicycles.size() - 1)
            tip_radius_out = fmaxf(e.a, .01);
    }

    return w(pos.real(), pos.imag(), z);
}

V estimate_normal(V p, J<V> &centers, V tip, F tip_radius) {
    I dummy;
    F eps = .001;
    I dummy_id, dummy_index;
    F dx = map(A(p, w(eps, 0, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(A(p, w(-eps, 0, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index);
    F dy = map(A(p, w(0, eps, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(A(p, w(0, -eps, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index);
    F dz = map(A(p, w(0, 0, eps)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(A(p, w(0, 0, -eps)), centers, tip, tip_radius, &dummy_id, &dummy_index);

    return N(w(dx, dy, dz));
}

void write_ppm_header(FILE *f) {
    fprintf(f, "P6\n%d %d\n255\n", W, H);
}

V get_light(F t) {
    F angle = t * 1.5;
    F radius = 3.;
    F x = cosf(angle) * radius;
    F y = sinf(angle * h) * 1.5;
    F z = sinf(angle) * radius;
    return N(w(x, y, z));
}

I main() {
    path = T(q, 64);
    compute_dft();
    for (I frame = 0; frame < FRAMES; ++frame) {
        char fname[64];
        snprintf(fname, sizeof(fname), "frame_%03d.ppm", frame);
        FILE *out = fopen(fname, "wb");
        if (!out) continue;
        write_ppm_header(out);

        F t = 2 * P * frame / FRAMES;
        F tip_radius;
        J<V> centers;
        V tip = get_tip_and_centers(t, centers, tip_radius);

        for (I y = 0; y < H; y++)
            for (I x = 0; x < W; x++)
                trail[y][x] *= .999;

        V light_dir = get_light(t);

        for (I y = 0; y < H; y++) {
            for (I x = 0; x < W; x++) {
                I shape_id = 0;
                I shape_index = -1;
    
                F u = (2. * x - W) / H;
                F v = (2. * y - H) / H;
                V ro = w(0, 0, -1);
                V rd = N(w(u, v, 1));
                F d = 0.;
                I id = 0;

                for (I i = 0; i < 100; i++) {
                    V p = A(ro, M(rd, d));
                    F dist = map(p, centers, tip, tip_radius, &shape_id, &shape_index);

                    if (dist < .001) { id = shape_id; break; }
                    if (d > 3.) break;
                    d += dist;
                }

                unsigned char r = 0, g = 0, b = 0;
                if (id == 1 || id == 2) {
                    V p = A(ro, M(rd, d));
                    V n = estimate_normal(p, centers, tip, tip_radius);
                    F ambient = .3;
                    F brightness = ambient + fmaxf(0., D(n, light_dir)) * .7;

                    F ball_t = (shape_index >= 0) ? (F)shape_index / (F)(centers.size() - 1) : 0.;

                    if (id == 1) {
                        F animated_t = fmodf(sqrt(ball_t) + frame * .01, 1.);
                        V color = U(animated_t);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                    } else if (id == 2) {
                        F tip_t = fmodf(sqrt(1.) + frame * .01, 1.);  // which is just fmodf(1.0 + frame * 0.01f, 1.0f)

                        V color = U(tip_t);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                        I tx = (int)((u + 1) * h * W);
                        I ty = (int)((v + 1) * h * H);
                        if (tx >= 0 && tx < W && ty >= 0 && ty < H)
                            trail[ty][tx] = 1.;
                    }
                }

                F trail_strength = trail[y][x];
                if (trail_strength > .001) {
                    // Mix original color with red trail
                    F tr = trail_strength * 255.;
                    r = (unsigned char)clamp((F)r + tr, 0.f, 255.f);
                    g = (unsigned char)clamp((F)g * (1. - trail_strength), 0., 255.);
                    b = (unsigned char)clamp((F)b * (1. - trail_strength), 0., 255.);
                }

                fputc(r, out); fputc(g, out); fputc(b, out);
            }
        }

        fclose(out);
    }

    return 0;
}