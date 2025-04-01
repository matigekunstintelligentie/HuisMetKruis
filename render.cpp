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

const F PI = 3.141592653589793f;

typedef struct {F x,y,z;} V;

V w(F x,F y,F z){return{x,y,z};}
V A(V a,V b){return w(a.x+b.x, a.y + b.y, a.z + b.z);}
V M(V a,F s){return w(a.x*s,a.y*s,a.z*s);}
V R(V a,V b){return w(a.x*b.x,a.y*b.y,a.z*b.z);}
F D(V a,V b){return a.x*b.x + a.y*b.y + a.z*b.z;}
F L(V v){return sqrtf(D(v,v));}
V N(V v){F l=L(v);return(l>0)?M(v,1.0f/l):v;}

struct Epicycle {
    I freq;
    F amp;
    F phase;
    complex<float> coeff;
};

V iq_palette(F t, V a, V b, V c, V d) {
    F x = cosf(2.0f * PI * (c.x * t + d.x));
    F y = cosf(2.0f * PI * (c.y * t + d.y));
    F z = cosf(2.0f * PI * (c.z * t + d.z));
    return A(a, R(b, w(x, y, z)));
}

V PALETTE_A = {0.5f, 0.5f, 0.5f};
V PALETTE_B = {0.5f, 0.5f, 0.5f};
V PALETTE_C = {1.0f, 1.0f, 1.0f};
V PALETTE_D = {0.0f, 0.33f, 0.67f};

vector<complex<float>> original_path = {
    {-0.5f, 0.5f},
    {-0.5f,  -0.5f},
    { 0.0f,  -0.85f},
    { 0.5f,  -0.5f},
    { 0.5f, 0.5f},
    {-0.5f, 0.5f},
    { 0.5f,  -0.5f},
    {-0.5f,  -0.5f},
    {0.5f, 0.5f},
    {-0.5f, 0.5f},
};

vector<complex<float>> interpolate_path(const vector<complex<float>>& points, I total_samples) {
    F total_length = 0.0f;
    vector<float> segment_lengths;
    for (size_t i = 0; i < points.size() - 1; ++i) {
        F len = abs(points[i + 1] - points[i]);
        segment_lengths.push_back(len);
        total_length += len;
    }

    vector<complex<float>> sampled_path;
    F step = total_length / total_samples;
    F current_dist = 0.0f;
    size_t seg = 0;
    F seg_pos = 0.0f;

    for (I i = 0; i < total_samples; ++i) {
        while (seg < segment_lengths.size() && seg_pos + segment_lengths[seg] < current_dist) {
            seg_pos += segment_lengths[seg];
            seg++;
        }
        if (seg >= segment_lengths.size()) break;
        F local_t = (current_dist - seg_pos) / segment_lengths[seg];
        complex<float> a = points[seg];
        complex<float> b = points[seg + 1];
        complex<float> interp = a + (b - a) * local_t;
        sampled_path.push_back(interp);
        current_dist += step;
    }

    return sampled_path;
}

vector<complex<float>> path;
vector<Epicycle> epicycles;

void compute_dft() {
    I N = path.size();
    epicycles.clear();
    for (I k = 0; k < N; ++k) {
        complex<float> sum(0.0f, 0.0f);
        for (I n = 0; n < N; ++n) {
            F phi = 2 * PI * k * n / N;
            sum += path[n] * exp(complex<float>(0, -phi));
        }
        sum /= static_cast<float>(N);
        I freq = (k <= N/2) ? k : k - N;
        epicycles.push_back({freq, abs(sum), atan2(sum.imag(), sum.real()), sum});
    }
    sort(epicycles.begin(), epicycles.end(), [](Epicycle a, Epicycle b) {
        return a.amp > b.amp;
    });
}

F trail[H][W] = {0};

F sdSphere(V p, V center, F r) {
    return L(A(p, {-center.x,-center.y,-center.z})) - r;
}

F map(V p, vector<V> &centers, V tip, F tip_radius, I *out_id, I *out_index) {
    F min_dist = 1e9f;
    I id = -1;
    I closest_index = -1;

    for (I i = 0; i < centers.size(); ++i) {
        F radius = fmaxf(epicycles[i].amp, 0.01f);
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

V get_tip_and_centers(F t, vector<V> &centers_out, F &tip_radius_out) {
    complex<float> pos(0.0f, 0.0f);
    centers_out.clear();

    F total_radius = 0.0f;
    for (const auto &e : epicycles)
        total_radius += e.amp;

    F z = total_radius + 0.1f;
    F dz = total_radius / epicycles.size();

    for (I i = 0; i < epicycles.size(); ++i) {
        Epicycle &e = epicycles[i];
        F angle = e.freq * t + e.phase;
        complex<float> offset = e.amp * exp(complex<float>(0, angle));
        centers_out.push_back(w(pos.real(), pos.imag(), z));
        pos += offset;
        z -= dz;
        if (i == epicycles.size() - 1)
            tip_radius_out = fmaxf(e.amp, 0.01f);
    }

    return w(pos.real(), pos.imag(), z);
}

V estimate_normal(V p, vector<V> &centers, V tip, F tip_radius) {
    I dummy;
    F eps = 0.001f;
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
    F angle = t * 1.5f;
    F radius = 3.0f;
    F x = cosf(angle) * radius;
    F y = sinf(angle * 0.5f) * 1.5f;
    F z = sinf(angle) * radius;
    return N(w(x, y, z));
}

I main() {
    path = interpolate_path(original_path, 64);
    compute_dft();
    for (I frame = 0; frame < FRAMES; ++frame) {
        char fname[64];
        snprintf(fname, sizeof(fname), "frame_%03d.ppm", frame);
        FILE *out = fopen(fname, "wb");
        if (!out) continue;
        write_ppm_header(out);

        F t = 2 * PI * frame / FRAMES;
        F tip_radius;
        vector<V> centers;
        V tip = get_tip_and_centers(t, centers, tip_radius);

        for (I y = 0; y < H; y++)
            for (I x = 0; x < W; x++)
                trail[y][x] *= 0.999f;

        V light_dir = get_light(t);

        for (I y = 0; y < H; y++) {
            for (I x = 0; x < W; x++) {
                I shape_id = 0;
                I shape_index = -1;
    
                F u = (2.0f * x - W) / H;
                F v = (2.0f * y - H) / H;
                V ro = w(0, 0, -1);
                V rd = N(w(u, v, 1));
                F d = 0.0f;
                I id = 0;

                for (I i = 0; i < 100; i++) {
                    V p = A(ro, M(rd, d));
                    F dist = map(p, centers, tip, tip_radius, &shape_id, &shape_index);

                    if (dist < 0.001f) { id = shape_id; break; }
                    if (d > 3.0f) break;
                    d += dist;
                }

                unsigned char r = 0, g = 0, b = 0;
                if (id == 1 || id == 2) {
                    V p = A(ro, M(rd, d));
                    V n = estimate_normal(p, centers, tip, tip_radius);
                    F ambient = 0.3f;
                    F brightness = ambient + fmaxf(0.0f, D(n, light_dir)) * 0.7f;

                    F ball_t = (shape_index >= 0) ? (float)shape_index / (float)(centers.size() - 1) : 0.0f;

                    if (id == 1) {
                        F animated_t = fmodf(sqrt(ball_t) + frame * 0.01f, 1.0f);
                        V color = iq_palette(animated_t, PALETTE_A, PALETTE_B, PALETTE_C, PALETTE_D);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                    } else if (id == 2) {
                        F tip_t = fmodf(sqrt(1.0f) + frame * 0.01f, 1.0f);  // which is just fmodf(1.0 + frame * 0.01f, 1.0f)

                        V color = iq_palette(tip_t, PALETTE_A, PALETTE_B, PALETTE_C, PALETTE_D);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                        I tx = (int)((u + 1) * 0.5f * W);
                        I ty = (int)((v + 1) * 0.5f * H);
                        if (tx >= 0 && tx < W && ty >= 0 && ty < H)
                            trail[ty][tx] = 1.0f;
                    }
                }

                F trail_strength = trail[y][x];
                if (trail_strength > 0.001f) {
                    // Mix original color with red trail
                    F tr = trail_strength * 255.0f;
                    r = (unsigned char)clamp((float)r + tr, 0.0f, 255.0f);
                    g = (unsigned char)clamp((float)g * (1.0f - trail_strength), 0.0f, 255.0f);
                    b = (unsigned char)clamp((float)b * (1.0f - trail_strength), 0.0f, 255.0f);
                }

                fputc(r, out); fputc(g, out); fputc(b, out);
            }
        }

        fclose(out);
    }

    return 0;
}