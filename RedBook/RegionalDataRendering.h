#pragma once

#include <glad/glad.h>
#include <glm.hpp>
#include <vector>

class RegionalDataRendering
{
public:
    void Init();
    void RenderBandedHeight(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap);
    void RenderHillshade(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap, float pixelWidthInMeters);
    void RenderRivers(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap, const float* depthmap, float pixelWidthInMeters);
    void Cleanup();
    bool IsInitialized() const { return initialzed; }
private:
    bool initialzed = false;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint bandedHeightProgram = 0;
    GLuint hillshadeProgram = 0;
    GLuint riversProgram = 0;
};

