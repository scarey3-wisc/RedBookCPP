#pragma once

#include <glad/glad.h>
#include <glm.hpp>
#include <vector>
class MeshPointRendering
{
public:
    void Init();
    void Render(std::vector<glm::vec2>& locs, float radius, float innerRadius, const glm::mat4& viewProj, const glm::vec4& color);
    void Cleanup();
    bool IsInitialized() const { return initialzed; }
private:
    bool initialzed = false;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint shaderProgram = 0;
    GLuint instanceVBO = 0;
};