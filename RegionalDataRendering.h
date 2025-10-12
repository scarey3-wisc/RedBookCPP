#pragma once

#include <glad/glad.h>
#include <glm.hpp>
#include <vector>

class RegionalDataRendering
{
public:
    void Init();
    void Render(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap, const glm::vec4& color);
    void Cleanup();
    bool IsInitialized() const { return initialzed; }
private:
    bool initialzed = false;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint shaderProgram = 0;
};

