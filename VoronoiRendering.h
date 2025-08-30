#pragma once

#include <glad/glad.h>
#include <glm.hpp>
#include <vector>
class VoronoiRendering
{
public:
    void Init();
    void Render(std::vector<glm::vec2>& locs, std::vector<glm::vec3>& colors, float range, const glm::mat4& viewProj);
    void Cleanup();
    bool IsInitialized() const { return initialzed; }
private:
    bool initialzed = false;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint shaderProgram = 0;
    GLuint instanceVBO = 0;
    GLuint colorVBO = 0;
};