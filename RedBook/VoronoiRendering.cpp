#include "VoronoiRendering.h"

#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include "RenderUtils.h"

static float rectVertices[] = {
    0.5f, -0.5f,
    -0.5f, -0.5f,
    0.5f, 0.5f,
    -0.5f, 0.5f
};

static const char* vertexSource = R"(
#version 330 core

layout (location = 0) in vec2 aPos;      // Rectangle vertex
layout (location = 1) in vec2 aOffset;   // Instance position
layout (location = 2) in vec3 aColor;    // Instance color

uniform mat4 uViewProj; // Your orthographic projection
uniform mat4 scale;

out vec3 fragLocalColor; // Pass color to fragment shader
out vec2 fragLocalPos;  // Pass to fragment shader

void main()
{
    fragLocalColor = aColor; // Pass the color to fragment shader
    fragLocalPos = aPos; // Pass the rectangle vertex to fragment shader
    gl_Position = uViewProj * (scale * vec4(aPos, 0.0, 1.0) + vec4(aOffset, 0.0, 0.0));
}
)";

static const char* fragmentSource = R"(
#version 330 core
in vec2 fragLocalPos; // Rectangle vertex from vertex shader
in vec3 fragLocalColor; // Color from vertex shader
out vec4 FragColor;

void main() 
{
    gl_FragDepth = length(fragLocalPos);
    
    // Set the color of the fragment
    FragColor = vec4(fragLocalColor, 1.0);
}
)";

void
VoronoiRendering::Init() {
    if (initialzed) {
        Cleanup(); // Clean up previous resources if already initialized
    }
    initialzed = true;

    // Create VAO and VBO
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);


    // Position attribute (layout location 0)
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(rectVertices), rectVertices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);


    // Offset attribute (layout location 1)
    glGenBuffers(1, &instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 10, nullptr, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glVertexAttribDivisor(1, 1); // Advance this attribute once per instance!


    // Color attribute (layout location 2)
    glGenBuffers(1, &colorVBO);
    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferData(GL_ARRAY_BUFFER, 10, nullptr, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glVertexAttribDivisor(2, 1); // Advance this attribute once per instance!

    glBindVertexArray(0);

    // Compile shaders (vertex & fragment)
    // You'll need to write helper functions to compile/link shaders
    shaderProgram = CompileAndLinkShader(vertexSource, fragmentSource);
}

void
VoronoiRendering::Render(std::vector<glm::vec2>& locs, std::vector<glm::vec3>& colors, float range, const glm::mat4& viewProj)
{

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glUseProgram(shaderProgram);

    // Build model matrix: translate and scale unit rectangle
    glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(range, range, 1.0f));


    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, locs.size() * sizeof(glm::vec2), locs.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), colors.data(), GL_DYNAMIC_DRAW);

    glm::mat4 mvp = viewProj;

    // Set uniforms
    GLint loc = glGetUniformLocation(shaderProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

    loc = glGetUniformLocation(shaderProgram, "scale");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));

    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, locs.size());
    glBindVertexArray(0);

    glDisable(GL_DEPTH_TEST);
}

void
VoronoiRendering::Cleanup() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (instanceVBO) glDeleteBuffers(1, &instanceVBO);
	if (colorVBO) glDeleteBuffers(1, &colorVBO);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (shaderProgram) glDeleteProgram(shaderProgram);
}