#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <chrono>
#include <algorithm>

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

constexpr int IX(int x, int y, int z, int N) {
    return x + y * N + z * N * N;
}

struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;
    int iter = 4;
    float buoyancy = 0.08f;
    float dissipation = 0.995f;
    float vorticity = 0.6f;
    float exposure = 0.4f;
    float maxVel = 0.5f;

    std::vector<float> s;
    std::vector<float> density;
    std::vector<float> curlX;
    std::vector<float> curlY;
    std::vector<float> curlZ;
    std::vector<float> curlM;

    std::vector<float> Vx;
    std::vector<float> Vy;
    std::vector<float> Vz;

    std::vector<float> Vx0;
    std::vector<float> Vy0;
    std::vector<float> Vz0;

    FluidCube(int size, float dt, float diff, float visc)
        : size(size), dt(dt), diff(diff), visc(visc),
          s(size * size * size, 0.0f),
          density(size * size * size, 0.0f),
          curlX(size * size * size, 0.0f),
          curlY(size * size * size, 0.0f),
          curlZ(size * size * size, 0.0f),
          curlM(size * size * size, 0.0f),
          Vx(size * size * size, 0.0f),
          Vy(size * size * size, 0.0f),
          Vz(size * size * size, 0.0f),
          Vx0(size * size * size, 0.0f),
          Vy0(size * size * size, 0.0f),
          Vz0(size * size * size, 0.0f)
    {}

    void addDensity(int x, int y, int z, float amount) {
        int N = size;
        if (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) return;
        density[IX(x, y, z, N)] += amount;
    }

    void addVelocity(int x, int y, int z, float ax, float ay, float az) {
        int N = size;
        if (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) return;
        int i = IX(x, y, z, N);
        Vx[i] += ax; Vy[i] += ay; Vz[i] += az;
    }

    void addDensitySphere(float cx, float cy, float cz, float r, float amount) {
        int x0 = (int)std::floor(cx - r), x1 = (int)std::ceil(cx + r);
        int y0 = (int)std::floor(cy - r), y1 = (int)std::ceil(cy + r);
        int z0 = (int)std::floor(cz - r), z1 = (int)std::ceil(cz + r);
        for (int z = z0; z <= z1; z++)
            for (int y = y0; y <= y1; y++)
                for (int x = x0; x <= x1; x++) {
                    float dx = x - cx, dy = y - cy, dz = z - cz;
                    float d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 <= r * r) addDensity(x, y, z, amount * (1.0f - d2 / (r * r + 1e-6f)));
                }
    }

    void addVelocitySphere(float cx, float cy, float cz, float r, float ax, float ay, float az) {
        int x0 = (int)std::floor(cx - r), x1 = (int)std::ceil(cx + r);
        int y0 = (int)std::floor(cy - r), y1 = (int)std::ceil(cy + r);
        int z0 = (int)std::floor(cz - r), z1 = (int)std::ceil(cz + r);
        for (int z = z0; z <= z1; z++)
            for (int y = y0; y <= y1; y++)
                for (int x = x0; x <= x1; x++) {
                    float dx = x - cx, dy = y - cy, dz = z - cz;
                    if (dx * dx + dy * dy + dz * dz <= r * r) addVelocity(x, y, z, ax, ay, az);
                }
    }

    void reset() {
        std::fill(density.begin(), density.end(), 0.0f);
        std::fill(Vx.begin(), Vx.end(), 0.0f);
        std::fill(Vy.begin(), Vy.end(), 0.0f);
        std::fill(Vz.begin(), Vz.end(), 0.0f);
    }

    void blendVelocitySphere(float cx, float cy, float cz, float r,
                             float ax, float ay, float az, float t) {
        int N = size;
        int x0 = (int)std::floor(cx - r), x1 = (int)std::ceil(cx + r);
        int y0 = (int)std::floor(cy - r), y1 = (int)std::ceil(cy + r);
        int z0 = (int)std::floor(cz - r), z1 = (int)std::ceil(cz + r);
        for (int z = z0; z <= z1; z++)
            for (int y = y0; y <= y1; y++)
                for (int x = x0; x <= x1; x++) {
                    if (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) continue;
                    float dx = x - cx, dy = y - cy, dz = z - cz;
                    if (dx * dx + dy * dy + dz * dz > r * r) continue;
                    int id = IX(x, y, z, N);
                    Vx[id] += (ax - Vx[id]) * t;
                    Vy[id] += (ay - Vy[id]) * t;
                    Vz[id] += (az - Vz[id]) * t;
                }
    }

    void step() {
        int N = size;
        vorticityConfinement();
        diffuse(1, Vx0.data(), Vx.data(), visc);
        diffuse(2, Vy0.data(), Vy.data(), visc);
        diffuse(3, Vz0.data(), Vz.data(), visc);

        project(Vx0.data(), Vy0.data(), Vz0.data(), Vx.data(), Vy.data());

        advect(1, Vx.data(), Vx0.data(), Vx0.data(), Vy0.data(), Vz0.data());
        advect(2, Vy.data(), Vy0.data(), Vx0.data(), Vy0.data(), Vz0.data());
        advect(3, Vz.data(), Vz0.data(), Vx0.data(), Vy0.data(), Vz0.data());

        project(Vx.data(), Vy.data(), Vz.data(), Vx0.data(), Vy0.data());

        diffuse(0, s.data(), density.data(), diff);
        advect(0, density.data(), s.data(), Vx.data(), Vy.data(), Vz.data());

        int NN = N * N;
        float buoy = buoyancy;
        float diss = dissipation;
        float mv = maxVel;
        float* d = density.data();
        float* vx = Vx.data();
        float* vy = Vy.data();
        float* vz = Vz.data();
        #pragma omp parallel for
        for (int i = 0; i < N * NN; i++) {
            vy[i] += buoy * d[i] * dt;
            d[i] *= diss;
            vx[i] = std::min(std::max(vx[i], -mv), mv);
            vy[i] = std::min(std::max(vy[i], -mv), mv);
            vz[i] = std::min(std::max(vz[i], -mv), mv);
        }
    }

private:
    void vorticityConfinement() {
        int N = size;
        float h = 0.5f;
        const float* u = Vx.data();
        const float* v = Vy.data();
        const float* w = Vz.data();
        int NN = N * N;
        #pragma omp parallel for
        for (int z = 1; z < N - 1; z++) {
            int zN = z * NN;
            for (int y = 1; y < N - 1; y++) {
                int base = zN + y * N;
                for (int i = 1; i < N - 1; i++) {
                    int id = base + i;
                    float cx = h * (w[id + N] - w[id - N] - v[id + NN] + v[id - NN]);
                    float cy = h * (u[id + NN] - u[id - NN] - w[id + 1] + w[id - 1]);
                    float cz = h * (v[id + 1] - v[id - 1] - u[id + N] + u[id - N]);
                    curlX[id] = cx;
                    curlY[id] = cy;
                    curlZ[id] = cz;
                    curlM[id] = std::sqrt(cx * cx + cy * cy + cz * cz);
                }
            }
        }
        float eps = vorticity;
        float* uf = Vx.data();
        float* vf = Vy.data();
        float* wf = Vz.data();
        const float* m = curlM.data();
        #pragma omp parallel for
        for (int z = 1; z < N - 1; z++) {
            int zN = z * NN;
            for (int y = 1; y < N - 1; y++) {
                int base = zN + y * N;
                for (int i = 1; i < N - 1; i++) {
                    int id = base + i;
                    float gx = h * (m[id + 1] - m[id - 1]);
                    float gy = h * (m[id + N] - m[id - N]);
                    float gz = h * (m[id + NN] - m[id - NN]);
                    float len = std::sqrt(gx * gx + gy * gy + gz * gz) + 1e-5f;
                    gx /= len; gy /= len; gz /= len;
                    float ox = curlX[id], oy = curlY[id], oz = curlZ[id];
                    uf[id] += eps * dt * (gy * oz - gz * oy);
                    vf[id] += eps * dt * (gz * ox - gx * oz);
                    wf[id] += eps * dt * (gx * oy - gy * ox);
                }
            }
        }
    }

    static void set_bnd(int b, float* x, int N) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(0, i, j, N)]     = b == 1 ? -x[IX(1, i, j, N)]     : x[IX(1, i, j, N)];
                x[IX(N - 1, i, j, N)] = b == 1 ? -x[IX(N - 2, i, j, N)] : x[IX(N - 2, i, j, N)];
                x[IX(i, 0, j, N)]     = b == 2 ? -x[IX(i, 1, j, N)]     : x[IX(i, 1, j, N)];
                x[IX(i, N - 1, j, N)] = b == 2 ? -x[IX(i, N - 2, j, N)] : x[IX(i, N - 2, j, N)];
                x[IX(i, j, 0, N)]     = b == 3 ? -x[IX(i, j, 1, N)]     : x[IX(i, j, 1, N)];
                x[IX(i, j, N - 1, N)] = b == 3 ? -x[IX(i, j, N - 2, N)] : x[IX(i, j, N - 2, N)];
            }
        }
        // edges: average of the two adjacent interior edge cells
        for (int k = 1; k < N - 1; k++) {
            x[IX(0, 0, k, N)]         = 0.5f * (x[IX(1, 0, k, N)]         + x[IX(0, 1, k, N)]);
            x[IX(N - 1, 0, k, N)]     = 0.5f * (x[IX(N - 2, 0, k, N)]     + x[IX(N - 1, 1, k, N)]);
            x[IX(0, N - 1, k, N)]     = 0.5f * (x[IX(1, N - 1, k, N)]     + x[IX(0, N - 2, k, N)]);
            x[IX(N - 1, N - 1, k, N)] = 0.5f * (x[IX(N - 2, N - 1, k, N)] + x[IX(N - 1, N - 2, k, N)]);
            x[IX(0, k, 0, N)]         = 0.5f * (x[IX(1, k, 0, N)]         + x[IX(0, k, 1, N)]);
            x[IX(N - 1, k, 0, N)]     = 0.5f * (x[IX(N - 2, k, 0, N)]     + x[IX(N - 1, k, 1, N)]);
            x[IX(0, k, N - 1, N)]     = 0.5f * (x[IX(1, k, N - 1, N)]     + x[IX(0, k, N - 2, N)]);
            x[IX(N - 1, k, N - 1, N)] = 0.5f * (x[IX(N - 2, k, N - 1, N)] + x[IX(N - 1, k, N - 2, N)]);
            x[IX(k, 0, 0, N)]         = 0.5f * (x[IX(k, 1, 0, N)]         + x[IX(k, 0, 1, N)]);
            x[IX(k, N - 1, 0, N)]     = 0.5f * (x[IX(k, N - 2, 0, N)]     + x[IX(k, N - 1, 1, N)]);
            x[IX(k, 0, N - 1, N)]     = 0.5f * (x[IX(k, 1, N - 1, N)]     + x[IX(k, 0, N - 2, N)]);
            x[IX(k, N - 1, N - 1, N)] = 0.5f * (x[IX(k, N - 2, N - 1, N)] + x[IX(k, N - 1, N - 2, N)]);
        }
        // corners: average of three adjacent cells
        x[IX(0, 0, 0, N)]             = (x[IX(1, 0, 0, N)]         + x[IX(0, 1, 0, N)]         + x[IX(0, 0, 1, N)])         / 3.0f;
        x[IX(N - 1, 0, 0, N)]         = (x[IX(N - 2, 0, 0, N)]     + x[IX(N - 1, 1, 0, N)]     + x[IX(N - 1, 0, 1, N)])     / 3.0f;
        x[IX(0, N - 1, 0, N)]         = (x[IX(1, N - 1, 0, N)]     + x[IX(0, N - 2, 0, N)]     + x[IX(0, N - 1, 1, N)])     / 3.0f;
        x[IX(N - 1, N - 1, 0, N)]     = (x[IX(N - 2, N - 1, 0, N)] + x[IX(N - 1, N - 2, 0, N)] + x[IX(N - 1, N - 1, 1, N)]) / 3.0f;
        x[IX(0, 0, N - 1, N)]         = (x[IX(1, 0, N - 1, N)]     + x[IX(0, 1, N - 1, N)]     + x[IX(0, 0, N - 2, N)])     / 3.0f;
        x[IX(N - 1, 0, N - 1, N)]     = (x[IX(N - 2, 0, N - 1, N)] + x[IX(N - 1, 1, N - 1, N)] + x[IX(N - 1, 0, N - 2, N)]) / 3.0f;
        x[IX(0, N - 1, N - 1, N)]     = (x[IX(1, N - 1, N - 1, N)] + x[IX(0, N - 2, N - 1, N)] + x[IX(0, N - 1, N - 2, N)]) / 3.0f;
        x[IX(N - 1, N - 1, N - 1, N)] = (x[IX(N - 2, N - 1, N - 1, N)] + x[IX(N - 1, N - 2, N - 1, N)] + x[IX(N - 1, N - 1, N - 2, N)]) / 3.0f;
    }

    void lin_solve(int b, float* x, const float* x0, float a, float c) {
        int N = size;
        float cInv = 1.0f / c;
        for (int k = 0; k < iter; k++) {
            for (int parity = 0; parity < 2; parity++) {
                #pragma omp parallel for
                for (int z = 1; z < N - 1; z++) {
                    int zN = z * N * N;
                    for (int y = 1; y < N - 1; y++) {
                        int base = zN + y * N;
                        for (int i = 1; i < N - 1; i++) {
                            if (((i + y + z) & 1) != parity) continue;
                            int id = base + i;
                            x[id] = (x0[id] + a * (x[id - 1] + x[id + 1] +
                                                   x[id - N] + x[id + N] +
                                                   x[id - N * N] + x[id + N * N])) * cInv;
                        }
                    }
                }
            }
            set_bnd(b, x, N);
        }
    }

    void diffuse(int b, float* x, const float* x0, float df) {
        float a = dt * df * (size - 2) * (size - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a);
    }

    void project(float* u, float* v, float* w, float* p, float* div) {
        int N = size;
        float h = 1.0f / (N - 2);
        #pragma omp parallel for
        for (int z = 1; z < N - 1; z++) {
            int zN = z * N * N;
            for (int y = 1; y < N - 1; y++) {
                int base = zN + y * N;
                for (int i = 1; i < N - 1; i++) {
                    int id = base + i;
                    div[id] = -0.5f * h * (u[id + 1] - u[id - 1] +
                                           v[id + N] - v[id - N] +
                                           w[id + N * N] - w[id - N * N]);
                    p[id] = 0.0f;
                }
            }
        }
        set_bnd(0, div, N);
        set_bnd(0, p, N);
        lin_solve(0, p, div, 1.0f, 6.0f);
        #pragma omp parallel for
        for (int z = 1; z < N - 1; z++) {
            int zN = z * N * N;
            for (int y = 1; y < N - 1; y++) {
                int base = zN + y * N;
                for (int i = 1; i < N - 1; i++) {
                    int id = base + i;
                    u[id] -= 0.5f * (p[id + 1] - p[id - 1]) / h;
                    v[id] -= 0.5f * (p[id + N] - p[id - N]) / h;
                    w[id] -= 0.5f * (p[id + N * N] - p[id - N * N]) / h;
                }
            }
        }
        set_bnd(1, u, N);
        set_bnd(2, v, N);
        set_bnd(3, w, N);
    }

    void advect(int b, float* d, const float* d0,
                const float* u, const float* v, const float* w) {
        int N = size;
        float dt0 = dt * (N - 2);
        #pragma omp parallel for
        for (int z = 1; z < N - 1; z++) {
            int zN = z * N * N;
            for (int y = 1; y < N - 1; y++) {
                int base = zN + y * N;
                for (int i = 1; i < N - 1; i++) {
                    int id = base + i;
                    float x = i - dt0 * u[id];
                    float yy = y - dt0 * v[id];
                    float zz = z - dt0 * w[id];
                    x = std::min(std::max(x, 0.5f), N - 1.5f);
                    yy = std::min(std::max(yy, 0.5f), N - 1.5f);
                    zz = std::min(std::max(zz, 0.5f), N - 1.5f);
                    int i0 = (int)x, i1 = i0 + 1;
                    int j0 = (int)yy, j1 = j0 + 1;
                    int k0 = (int)zz, k1 = k0 + 1;
                    float s1 = x - i0, s0 = 1.0f - s1;
                    float t1 = yy - j0, t0 = 1.0f - t1;
                    float u1 = zz - k0, u0 = 1.0f - u1;
                    int NN = N * N;
                    d[id] =
                        u0 * (t0 * (s0 * d0[i0 + j0 * N + k0 * NN] + s1 * d0[i1 + j0 * N + k0 * NN]) +
                              t1 * (s0 * d0[i0 + j1 * N + k0 * NN] + s1 * d0[i1 + j1 * N + k0 * NN])) +
                        u1 * (t0 * (s0 * d0[i0 + j0 * N + k1 * NN] + s1 * d0[i1 + j0 * N + k1 * NN]) +
                              t1 * (s0 * d0[i0 + j1 * N + k1 * NN] + s1 * d0[i1 + j1 * N + k1 * NN]));
                }
            }
        }
        set_bnd(b, d, N);
    }
};

static const char* volumeVS = R"GLSL(
#version 330

in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;

uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;

out vec3 fragPosition;
out vec2 fragTexCoord;
out vec3 fragNormal;
out vec4 fragColor;

void main() {
    vec4 worldPos = matModel * vec4(vertexPosition, 1.0);
    fragPosition = worldPos.xyz;
    fragTexCoord = vertexTexCoord;
    fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 0.0)));
    fragColor = vertexColor;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
)GLSL";

static const char* volumeFS = R"GLSL(
#version 330

in vec3 fragPosition;
in vec2 fragTexCoord;
in vec3 fragNormal;
in vec4 fragColor;

out vec4 finalColor;

uniform sampler2D atlas;
uniform vec3 camPos;
uniform vec3 boxMin;
uniform vec3 boxMax;
uniform int tilesX;
uniform int tilesY;
uniform int gridN;
uniform int colorMode;
uniform float marchSteps;

float sampleDensity(vec3 p) {
    p = clamp(p, 0.0, 0.99999);
    float zf = clamp(p.z * float(gridN) - 0.5, 0.0, float(gridN - 1));
    int z0 = int(floor(zf));
    int z1 = min(z0 + 1, gridN - 1);
    float fz = zf - float(z0);
    vec2 inSlice = clamp(p.xy * float(gridN), 0.5, float(gridN) - 0.5);
    float N = float(gridN);
    vec2 atlasSize = vec2(float(tilesX) * N, float(tilesY) * N);
    vec2 uv0 = (vec2(float(z0 - (z0 / tilesX) * tilesX), float(z0 / tilesX)) * N + inSlice) / atlasSize;
    vec2 uv1 = (vec2(float(z1 - (z1 / tilesX) * tilesX), float(z1 / tilesX)) * N + inSlice) / atlasSize;
    float a = texture(atlas, uv0).r;
    float b = texture(atlas, uv1).r;
    return mix(a, b, fz);
}

vec3 ramp(float d, float acc) {
    if (colorMode == 0) {
        return mix(vec3(0.55, 0.6, 0.7), vec3(1.0), clamp(d * 1.5, 0.0, 1.0));
    } else if (colorMode == 1) {
        vec3 c = mix(vec3(0.1, 0.15, 0.5), vec3(1.0, 0.45, 0.05), clamp(d * 2.0, 0.0, 1.0));
        c = mix(c, vec3(1.0, 0.95, 0.8), clamp(d * 2.0 - 1.0, 0.0, 1.0));
        return c * 1.6;
    } else {
        float h = fract(acc * 0.7 + d);
        return clamp(abs(mod(h * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    }
}

vec2 boxIntersect(vec3 ro, vec3 rd) {
    vec3 inv = 1.0 / rd;
    vec3 t0 = (boxMin - ro) * inv;
    vec3 t1 = (boxMax - ro) * inv;
    vec3 tmin = min(t0, t1), tmax = max(t0, t1);
    return vec2(max(max(tmin.x, tmin.y), tmin.z), min(min(tmax.x, tmax.y), tmax.z));
}

void main() {
    vec3 rd = normalize(fragPosition - camPos);
    vec2 t = boxIntersect(camPos, rd);
    float t0 = max(t.x, 0.0);
    float t1 = t.y;
    if (t1 <= t0) discard;

    float steps = marchSteps;
    float stepLen = (t1 - t0) / steps;
    float sigma = 90.0;
    float T = 1.0;
    vec3 col = vec3(0.0);
    vec3 size = boxMax - boxMin;

    for (int i = 0; i < 256; i++) {
        if (float(i) >= steps) break;
        vec3 pos = camPos + rd * (t0 + (float(i) + 0.5) * stepLen);
        vec3 uvw = (pos - boxMin) / size;
        float d = sampleDensity(uvw);
        if (d > 0.003) {
            float a = 1.0 - exp(-d * sigma * stepLen);
            col += T * a * ramp(d, 1.0 - T);
            T *= exp(-d * sigma * stepLen);
            if (T < 0.01) break;
        }
    }
    float alpha = 1.0 - T;
    if (alpha < 0.005) discard;
    finalColor = vec4(col, alpha);
}
)GLSL";

static Color sliceColor(float d, int mode) {
    d = std::min(std::max(d, 0.0f), 1.0f);
    if (mode == 0) {
        unsigned char v = (unsigned char)(d * 255);
        return Color{v, v, v, 255};
    }
    if (mode == 1) {
        float t = std::min(d * 2.0f, 1.0f);
        Color a{13, 26, 102, 255}, b{255, 115, 13, 255}, c{255, 242, 204, 255};
        Color lo = {(unsigned char)(a.r + (b.r - a.r) * t), (unsigned char)(a.g + (b.g - a.g) * t),
                    (unsigned char)(a.b + (b.b - a.b) * t), 255};
        float u = std::min(std::max(d * 2.0f - 1.0f, 0.0f), 1.0f);
        return Color{(unsigned char)(lo.r + (c.r - lo.r) * u), (unsigned char)(lo.g + (c.g - lo.g) * u),
                     (unsigned char)(lo.b + (c.b - lo.b) * u), 255};
    }
    return Color{(unsigned char)(d * 255), (unsigned char)(128 + 127 * std::sin(d * 6.28f)),
                 (unsigned char)(255 - d * 255), 255};
}

int main(int argc, char** argv) {
    bool shotVolume = false, shotSlice = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--screenshot") == 0) shotVolume = true;
        if (std::strcmp(argv[i], "--screenshot-slice") == 0) shotSlice = true;
    }

    const int N = 48;
    FluidCube fluid(N, 0.1f, 0.00001f, 0.00001f);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
    InitWindow(1280, 720, "3D Fluid - Stable Fluids (raylib)");
    SetTargetFPS(60);

    int tilesX = (int)std::ceil(std::sqrt((double)N));
    int tilesY = (N + tilesX - 1) / tilesX;
    int atlasW = tilesX * N, atlasH = tilesY * N;
    std::vector<unsigned char> atlasBuf(atlasW * atlasH, 0);

    Image atlasImg;
    atlasImg.data = atlasBuf.data();
    atlasImg.width = atlasW;
    atlasImg.height = atlasH;
    atlasImg.mipmaps = 1;
    atlasImg.format = PIXELFORMAT_UNCOMPRESSED_GRAYSCALE;
    Texture2D atlasTex = LoadTextureFromImage(atlasImg);
    SetTextureFilter(atlasTex, TEXTURE_FILTER_BILINEAR);
    atlasImg.data = nullptr;

    std::vector<Color> sliceBuf(N * N);
    Image sliceImg{sliceBuf.data(), N, N, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8};
    Texture2D sliceTex = LoadTextureFromImage(sliceImg);
    SetTextureFilter(sliceTex, TEXTURE_FILTER_BILINEAR);

    Shader volShader = LoadShaderFromMemory(volumeVS, volumeFS);
    bool shaderOk = volShader.id > 0 && volShader.id != rlGetShaderIdDefault();
    if (!shaderOk) TraceLog(LOG_WARNING, "Volume shader failed, falling back to slice view");

    int locAtlas = GetShaderLocation(volShader, "atlas");
    int locCamPos = GetShaderLocation(volShader, "camPos");
    int locBoxMin = GetShaderLocation(volShader, "boxMin");
    int locBoxMax = GetShaderLocation(volShader, "boxMax");
    int locTilesX = GetShaderLocation(volShader, "tilesX");
    int locTilesY = GetShaderLocation(volShader, "tilesY");
    int locGridN = GetShaderLocation(volShader, "gridN");
    int locColorMode = GetShaderLocation(volShader, "colorMode");
    int locMarch = GetShaderLocation(volShader, "marchSteps");

    Mesh cubeMesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model cubeModel = LoadModelFromMesh(cubeMesh);
    cubeModel.materials[0].shader = volShader;

    Vector3 boxMin{0, 0, 0}, boxMax{1, 1, 1};
    Vector3 boxCenter{0.5f, 0.5f, 0.5f};

    Camera3D cam{};
    cam.up = {0, 1, 0};
    cam.fovy = 50.0f;
    cam.projection = CAMERA_PERSPECTIVE;
    float yaw = 0.8f, pitch = 0.35f, dist = 2.6f;

    bool paused = false, emitter = true, sliceView = !shaderOk || shotSlice;
    int colorMode = 1;
    bool showHelp = true;
    float stepMs = 0.0f;
    double stepMsAccum = 0.0;
    int stepCount = 0;
    int frame = 0;
    bool mouseDown = false;
    Vector3 prevHit{};
    float emitAngle = 0.0f;

    auto updateCam = [&]() {
        float cp = std::cos(pitch), sp = std::sin(pitch);
        cam.position = {boxCenter.x + dist * cp * std::cos(yaw),
                        boxCenter.y + dist * sp,
                        boxCenter.z + dist * cp * std::sin(yaw)};
        cam.target = boxCenter;
    };

    auto mouseHitOnPlane = [&](Vector2 m, Vector3& out) -> bool {
        Ray ray = GetScreenToWorldRay(m, cam);
        Vector3 n = Vector3Normalize(Vector3Subtract(cam.target, cam.position));
        float denom = Vector3DotProduct(ray.direction, n);
        if (std::fabs(denom) < 1e-5f) return false;
        float t = Vector3DotProduct(Vector3Subtract(boxCenter, ray.position), n) / denom;
        if (t <= 0) return false;
        out = Vector3Add(ray.position, Vector3Scale(ray.direction, t));
        return true;
    };

    auto worldToGrid = [](Vector3 w) {
        return Vector3{w.x * N, w.y * N, w.z * N};
    };

    while (!WindowShouldClose()) {
        frame++;
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 d = GetMouseDelta();
            yaw -= d.x * 0.005f;
            pitch += d.y * 0.005f;
            pitch = std::min(std::max(pitch, -1.5f), 1.5f);
        }
        dist -= GetMouseWheelMove() * 0.15f;
        dist = std::min(std::max(dist, 0.4f), 12.0f);
        updateCam();

        if (IsKeyPressed(KEY_P)) paused = !paused;
        if (IsKeyPressed(KEY_R)) fluid.reset();
        if (IsKeyPressed(KEY_C)) colorMode = (colorMode + 1) % 3;
        if (IsKeyPressed(KEY_V)) sliceView = !sliceView;
        if (IsKeyPressed(KEY_SPACE)) emitter = !emitter;
        if (IsKeyPressed(KEY_F1)) showHelp = !showHelp;
        if (IsKeyPressed(KEY_ONE)) fluid.buoyancy = std::max(0.0f, fluid.buoyancy - 0.05f);
        if (IsKeyPressed(KEY_TWO)) fluid.buoyancy += 0.05f;
        if (IsKeyPressed(KEY_THREE)) fluid.dissipation = std::max(0.9f, fluid.dissipation - 0.002f);
        if (IsKeyPressed(KEY_FOUR)) fluid.dissipation = std::min(1.0f, fluid.dissipation + 0.002f);
        if (IsKeyPressed(KEY_LEFT_BRACKET)) fluid.iter = std::max(1, fluid.iter - 1);
        if (IsKeyPressed(KEY_RIGHT_BRACKET)) fluid.iter = std::min(64, fluid.iter + 1);
        if (IsKeyPressed(KEY_FIVE)) fluid.vorticity = std::max(0.0f, fluid.vorticity - 0.1f);
        if (IsKeyPressed(KEY_SIX)) fluid.vorticity = std::min(10.0f, fluid.vorticity + 0.1f);
        if (IsKeyPressed(KEY_SEVEN)) fluid.exposure = std::max(0.01f, fluid.exposure - 0.05f);
        if (IsKeyPressed(KEY_EIGHT)) fluid.exposure = std::min(5.0f, fluid.exposure + 0.05f);

        if (!paused) {
            emitAngle += 0.05f;
            if (emitter) {
                float cx = N * 0.5f + std::cos(emitAngle) * 2.0f;
                float cz = N * 0.5f + std::sin(emitAngle) * 2.0f;
                fluid.addDensitySphere(cx, 5.0f, cz, 2.5f, 5.0f);
                fluid.blendVelocitySphere(cx, 5.0f, cz, 3.0f,
                                          std::cos(emitAngle * 1.7f) * 0.12f,
                                          0.28f,
                                          std::sin(emitAngle * 1.7f) * 0.12f, 0.5f);
            }

            Vector3 hit;
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && mouseHitOnPlane(GetMousePosition(), hit)) {
                Vector3 g = worldToGrid(hit);
                fluid.addDensitySphere(g.x, g.y, g.z, 3.0f, 8.0f);
                if (mouseDown) {
                    Vector3 drag = Vector3Subtract(hit, prevHit);
                    float scale = 12.0f;
                    fluid.addVelocitySphere(g.x, g.y, g.z, 4.0f,
                                            drag.x * scale, drag.y * scale, drag.z * scale);
                }
                prevHit = hit;
                mouseDown = true;
            } else {
                mouseDown = false;
            }

            auto t0 = std::chrono::high_resolution_clock::now();
            fluid.step();
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            stepMsAccum += ms;
            stepCount++;
            stepMs = (float)(stepMsAccum / stepCount);
        }

        float expo = fluid.exposure;
        for (int z = 0; z < N; z++) {
            for (int y = 0; y < N; y++) {
                int ty = z / tilesX, tx = z - ty * tilesX;
                unsigned char* row = atlasBuf.data() + (ty * N + y) * atlasW + tx * N;
                const float* src = fluid.density.data() + z * N * N + y * N;
                for (int x = 0; x < N; x++) {
                    float v = 1.0f - std::exp(-src[x] * expo);
                    row[x] = (unsigned char)(std::min(v, 1.0f) * 255.0f);
                }
            }
        }
        UpdateTexture(atlasTex, atlasBuf.data());

        int midZ = N / 2;
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++) {
                float v = 1.0f - std::exp(-fluid.density[IX(x, N - 1 - y, midZ, N)] * expo);
                sliceBuf[y * N + x] = sliceColor(v, colorMode);
            }
        UpdateTexture(sliceTex, sliceBuf.data());

        BeginDrawing();
        ClearBackground(Color{18, 20, 26, 255});

        if (!sliceView) {
            BeginMode3D(cam);
            DrawGrid(20, 0.25f);
            DrawCubeWiresV(boxCenter, {1, 1, 1}, Color{120, 140, 180, 255});
            if (shaderOk) {
                SetShaderValue(volShader, locCamPos, &cam.position, SHADER_UNIFORM_VEC3);
                SetShaderValue(volShader, locBoxMin, &boxMin, SHADER_UNIFORM_VEC3);
                SetShaderValue(volShader, locBoxMax, &boxMax, SHADER_UNIFORM_VEC3);
                SetShaderValue(volShader, locTilesX, &tilesX, SHADER_UNIFORM_INT);
                SetShaderValue(volShader, locTilesY, &tilesY, SHADER_UNIFORM_INT);
                SetShaderValue(volShader, locGridN, &N, SHADER_UNIFORM_INT);
                SetShaderValue(volShader, locColorMode, &colorMode, SHADER_UNIFORM_INT);
                float steps = 96.0f;
                SetShaderValue(volShader, locMarch, &steps, SHADER_UNIFORM_FLOAT);
                SetShaderValueTexture(volShader, locAtlas, atlasTex);
                rlEnableBackfaceCulling();
                rlSetCullFace(RL_CULL_FACE_FRONT);
                rlDisableDepthMask();
                DrawModel(cubeModel, boxCenter, 1.0f, WHITE);
                rlEnableDepthMask();
                rlSetCullFace(RL_CULL_FACE_BACK);
            }
            EndMode3D();
        } else {
            float top = 180.0f;
            float avail = GetScreenHeight() - top - 40.0f;
            float scale = std::min(GetScreenWidth() - 40.0f, avail) / (float)N;
            float ox = (GetScreenWidth() - N * scale) * 0.5f;
            float oy = top + (avail - N * scale) * 0.5f;
            DrawTextureEx(sliceTex, {ox, oy}, 0.0f, scale, WHITE);
            int mid = N / 2;
            for (int y = 0; y < N; y += 4)
                for (int x = 0; x < N; x += 4) {
                    int id = IX(x, N - 1 - y, mid, N);
                    float px = ox + (x + 0.5f) * scale, py = oy + (y + 0.5f) * scale;
                    DrawLine((int)px, (int)py,
                             (int)(px + fluid.Vx[id] * 80.0f), (int)(py - fluid.Vy[id] * 80.0f),
                             Color{140, 200, 255, 160});
                }
            DrawRectangleLines((int)ox, (int)oy, (int)(N * scale), (int)(N * scale), GRAY);
            DrawText("Slice view (V toggles back)", (int)ox, (int)(oy + N * scale + 8), 18, LIGHTGRAY);
        }

        DrawFPS(10, 10);
        DrawText(TextFormat("N=%d  step=%.2f ms  buoy=%.2f  diss=%.3f  iter=%d  vort=%.1f  expo=%.2f  mode=%d%s%s",
                            N, stepMs, fluid.buoyancy, fluid.dissipation, fluid.iter, fluid.vorticity,
                            fluid.exposure, colorMode,
                            paused ? "  [PAUSED]" : "", emitter ? "  [EMIT]" : ""),
                 10, 34, 16, RAYWHITE);
        if (showHelp) {
            int y = 56;
            const char* lines[] = {
                "LMB drag: inject smoke + velocity",
                "RMB drag: orbit   Wheel: zoom",
                "Space: emitter  P: pause  R: reset",
                "C: color mode  V: slice/volume view",
                "1/2: buoyancy  3/4: dissipation  [/]: solver iters",
                "5/6: vorticity  7/8: exposure",
                "F1: toggle help"
            };
            for (auto l : lines) { DrawText(l, 10, y, 14, Color{200, 205, 215, 255}); y += 17; }
        }
        EndDrawing();

        if (shotVolume && frame >= 180) {
            TakeScreenshot("shot.png");
            printf("avg step time: %.3f ms over %d steps\n", stepMsAccum / stepCount, stepCount);
            break;
        }
        if (shotSlice && frame >= 180) {
            TakeScreenshot("shot_slice.png");
            printf("avg step time: %.3f ms over %d steps\n", stepMsAccum / stepCount, stepCount);
            break;
        }
    }

    UnloadTexture(atlasTex);
    UnloadTexture(sliceTex);
    UnloadModel(cubeModel);
    UnloadShader(volShader);
    CloseWindow();
    return 0;
}
