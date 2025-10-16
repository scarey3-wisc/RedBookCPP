#pragma once

#include <string>

static const std::string sourceHeader = R"(
#version 330 core
)";

static const std::string hsb2rgb = R"(
vec3 hsb2rgb(in vec3 c) {
    vec3 rgb = clamp(
        abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0,
        0.0,
        1.0
    );
    rgb = rgb * rgb * (3.0 - 2.0 * rgb); // smooth interpolation
    return c.z * mix(vec3(1.0), rgb, c.y);
})";

static const std::string getHeightmapNormal = R"(
vec3 getHeightmapNormal(isampler2D heightTex, vec2 uv, float mpp, float int2float) {
    vec2 texSize = textureSize(heightTex, 0);
    float texelSize = 1.0 / texSize.x; // assuming square texture
    float hL = float(texture(heightTex, uv + vec2(-texelSize, 0.0)).r) / int2float;
    float hR = float(texture(heightTex, uv + vec2(texelSize, 0.0)).r) / int2float;
    float hD = float(texture(heightTex, uv + vec2(0.0, -texelSize)).r) / int2float;
    float hU = float(texture(heightTex, uv + vec2(0.0, texelSize)).r) / int2float;
    float dx = (hR - hL) * 0.5 / mpp;
    float dy = (hU - hD) * 0.5 / mpp;
    return normalize(vec3(-dx, -dy, 1.0));
}
)";