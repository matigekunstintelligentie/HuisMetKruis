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
const float PI = 3.141592653589793f;

typedef struct { float x, y, z; } Vec3;

Vec3 vec3(float x, float y, float z) { return {x, y, z}; }
Vec3 add(Vec3 a, Vec3 b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
Vec3 sub(Vec3 a, Vec3 b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
Vec3 mul(Vec3 a, float s) { return vec3(a.x * s, a.y * s, a.z * s); }

Vec3 mul_vec(Vec3 a, Vec3 b) {
    return vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}
float dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
float length3(Vec3 v) { return sqrtf(dot(v, v)); }
Vec3 normalize(Vec3 v) { float l = length3(v); return (l > 0) ? mul(v, 1.0f/l) : v; }

struct Epicycle {
    int freq;
    float amp;
    float phase;
    complex<float> coeff;
};

Vec3 iq_palette(float t, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    float x = cosf(2.0f * PI * (c.x * t + d.x));
    float y = cosf(2.0f * PI * (c.y * t + d.y));
    float z = cosf(2.0f * PI * (c.z * t + d.z));
    return add(a, mul_vec(b, vec3(x, y, z)));
}

Vec3 PALETTE_A = {0.5f, 0.5f, 0.5f};
Vec3 PALETTE_B = {0.5f, 0.5f, 0.5f};
Vec3 PALETTE_C = {1.0f, 1.0f, 1.0f};
Vec3 PALETTE_D = {0.0f, 0.33f, 0.67f};

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

vector<complex<float>> interpolate_path(const vector<complex<float>>& points, int total_samples) {
    float total_length = 0.0f;
    vector<float> segment_lengths;
    for (size_t i = 0; i < points.size() - 1; ++i) {
        float len = abs(points[i + 1] - points[i]);
        segment_lengths.push_back(len);
        total_length += len;
    }

    vector<complex<float>> sampled_path;
    float step = total_length / total_samples;
    float current_dist = 0.0f;
    size_t seg = 0;
    float seg_pos = 0.0f;

    for (int i = 0; i < total_samples; ++i) {
        while (seg < segment_lengths.size() && seg_pos + segment_lengths[seg] < current_dist) {
            seg_pos += segment_lengths[seg];
            seg++;
        }
        if (seg >= segment_lengths.size()) break;
        float local_t = (current_dist - seg_pos) / segment_lengths[seg];
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
    int N = path.size();
    epicycles.clear();
    for (int k = 0; k < N; ++k) {
        complex<float> sum(0.0f, 0.0f);
        for (int n = 0; n < N; ++n) {
            float phi = 2 * PI * k * n / N;
            sum += path[n] * exp(complex<float>(0, -phi));
        }
        sum /= static_cast<float>(N);
        int freq = (k <= N/2) ? k : k - N;
        epicycles.push_back({freq, abs(sum), atan2(sum.imag(), sum.real()), sum});
    }
    sort(epicycles.begin(), epicycles.end(), [](Epicycle a, Epicycle b) {
        return a.amp > b.amp;
    });
}

float trail[H][W] = {0};

float sdSphere(Vec3 p, Vec3 center, float r) {
    return length3(sub(p, center)) - r;
}

float map(Vec3 p, vector<Vec3> &centers, Vec3 tip, float tip_radius, int *out_id, int *out_index) {
    float min_dist = 1e9f;
    int id = -1;
    int closest_index = -1;

    for (int i = 0; i < centers.size(); ++i) {
        float radius = fmaxf(epicycles[i].amp, 0.01f);
        float d = sdSphere(p, centers[i], radius);
        if (d < min_dist) {
            min_dist = d;
            id = 1;
            closest_index = i;
        }
    }
    float d_tip = sdSphere(p, tip, tip_radius);
    if (d_tip < min_dist) {
        min_dist = d_tip;
        id = 2;
        closest_index = -1;
    }
    *out_id = id;
    *out_index = closest_index;
    return min_dist;
}

Vec3 get_tip_and_centers(float t, vector<Vec3> &centers_out, float &tip_radius_out) {
    complex<float> pos(0.0f, 0.0f);
    centers_out.clear();

    float total_radius = 0.0f;
    for (const auto &e : epicycles)
        total_radius += e.amp;

    float z = total_radius + 0.1f;
    float dz = total_radius / epicycles.size();

    for (int i = 0; i < epicycles.size(); ++i) {
        Epicycle &e = epicycles[i];
        float angle = e.freq * t + e.phase;
        complex<float> offset = e.amp * exp(complex<float>(0, angle));
        centers_out.push_back(vec3(pos.real(), pos.imag(), z));
        pos += offset;
        z -= dz;
        if (i == epicycles.size() - 1)
            tip_radius_out = fmaxf(e.amp, 0.01f);
    }

    return vec3(pos.real(), pos.imag(), z);
}

Vec3 estimate_normal(Vec3 p, vector<Vec3> &centers, Vec3 tip, float tip_radius) {
    int dummy;
    float eps = 0.001f;
    int dummy_id, dummy_index;
    float dx = map(add(p, vec3(eps, 0, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(add(p, vec3(-eps, 0, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index);
    float dy = map(add(p, vec3(0, eps, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(add(p, vec3(0, -eps, 0)), centers, tip, tip_radius, &dummy_id, &dummy_index);
    float dz = map(add(p, vec3(0, 0, eps)), centers, tip, tip_radius, &dummy_id, &dummy_index) -
               map(add(p, vec3(0, 0, -eps)), centers, tip, tip_radius, &dummy_id, &dummy_index);

    return normalize(vec3(dx, dy, dz));
}

void write_ppm_header(FILE *f) {
    fprintf(f, "P6\n%d %d\n255\n", W, H);
}

Vec3 get_light(float t) {
    float angle = t * 1.5f;
    float radius = 3.0f;
    float x = cosf(angle) * radius;
    float y = sinf(angle * 0.5f) * 1.5f;
    float z = sinf(angle) * radius;
    return normalize(vec3(x, y, z));
}

int main() {
    path = interpolate_path(original_path, 64);
    compute_dft();
    for (int frame = 0; frame < FRAMES; ++frame) {
        char fname[64];
        snprintf(fname, sizeof(fname), "frame_%03d.ppm", frame);
        FILE *out = fopen(fname, "wb");
        if (!out) continue;
        write_ppm_header(out);

        float t = 2 * PI * frame / FRAMES;
        float tip_radius;
        vector<Vec3> centers;
        Vec3 tip = get_tip_and_centers(t, centers, tip_radius);

        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++)
                trail[y][x] *= 0.999f;

        Vec3 light_dir = get_light(t);

        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                int shape_id = 0;
                int shape_index = -1;
    
                float u = (2.0f * x - W) / H;
                float v = (2.0f * y - H) / H;
                Vec3 ro = vec3(0, 0, -1);
                Vec3 rd = normalize(vec3(u, v, 1));
                float d = 0.0f;
                int id = 0;

                for (int i = 0; i < 100; i++) {
                    Vec3 p = add(ro, mul(rd, d));
                    float dist = map(p, centers, tip, tip_radius, &shape_id, &shape_index);

                    if (dist < 0.001f) { id = shape_id; break; }
                    if (d > 3.0f) break;
                    d += dist;
                }

                unsigned char r = 0, g = 0, b = 0;
                if (id == 1 || id == 2) {
                    Vec3 p = add(ro, mul(rd, d));
                    Vec3 n = estimate_normal(p, centers, tip, tip_radius);
                    float ambient = 0.3f;
                    float brightness = ambient + fmaxf(0.0f, dot(n, light_dir)) * 0.7f;

                    float ball_t = (shape_index >= 0) ? (float)shape_index / (float)(centers.size() - 1) : 0.0f;

                    if (id == 1) {
                        float animated_t = fmodf(sqrt(ball_t) + frame * 0.01f, 1.0f);
                        Vec3 color = iq_palette(animated_t, PALETTE_A, PALETTE_B, PALETTE_C, PALETTE_D);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                    } else if (id == 2) {
                        float tip_t = fmodf(sqrt(1.0f) + frame * 0.01f, 1.0f);  // which is just fmodf(1.0 + frame * 0.01f, 1.0f)

                        Vec3 color = iq_palette(tip_t, PALETTE_A, PALETTE_B, PALETTE_C, PALETTE_D);
                        r = (unsigned char)(brightness * 255 * color.x);
                        g = (unsigned char)(brightness * 255 * color.y);
                        b = (unsigned char)(brightness * 255 * color.z);
                        int tx = (int)((u + 1) * 0.5f * W);
                        int ty = (int)((v + 1) * 0.5f * H);
                        if (tx >= 0 && tx < W && ty >= 0 && ty < H)
                            trail[ty][tx] = 1.0f;
                    }
                }

                float trail_strength = trail[y][x];
                if (trail_strength > 0.001f) {
                    // Mix original color with red trail
                    float tr = trail_strength * 255.0f;
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