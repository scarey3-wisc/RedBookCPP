#include "MeshPointRendering.h"

#include "OutlineRendering.h"
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

uniform mat4 uViewProj; // Your orthographic projection
uniform mat4 scale;

out vec2 fragLocalPos;  // Pass to fragment shader

void main()
{
    fragLocalPos = aPos; // Pass the rectangle vertex to fragment shader
    gl_Position = uViewProj * (scale * vec4(aPos, 0.0, 1.0) + vec4(aOffset, 0.0, 0.0));
}
)";

static const char* fragmentSource = R"(
#version 330 core
in vec2 fragLocalPos; // Rectangle vertex from vertex shader
out vec4 FragColor;

uniform vec4 uColor;
uniform float inRad;

void main() 
{
    float dist = length(fragLocalPos);
    if (dist > 0.5)
        discard;  // outside circle
    if(dist < inRad)
        discard; // inside anulus
    // Set the color of the fragment
    FragColor = uColor;
}
)";

void
MeshPointRendering::Init() {
    if (initialzed) {
        Cleanup(); // Clean up previous resources if already initialized
    }
    initialzed = true;


    // Create VAO and VBO
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(rectVertices), rectVertices, GL_STATIC_DRAW);

    // Position attribute (layout location 0)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);

    glGenBuffers(1, &instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 10, nullptr, GL_DYNAMIC_DRAW);

    // Offset attribute (layout location 1)
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glVertexAttribDivisor(1, 1); // Advance this attribute once per instance!

    glBindVertexArray(0);

    // Compile shaders (vertex & fragment)
    // You'll need to write helper functions to compile/link shaders
    shaderProgram = CompileAndLinkShader(vertexSource, fragmentSource);
}

void
MeshPointRendering::Render(std::vector<glm::vec2>& locs, float radius, float innerRadius, const glm::mat4& viewProj, const glm::vec4& color)
{

    glUseProgram(shaderProgram);

    // Build model matrix: translate and scale unit rectangle
    glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(radius, radius, 1.0f));


    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, locs.size() * sizeof(glm::vec2), locs.data(), GL_DYNAMIC_DRAW);

    glm::mat4 mvp = viewProj;

	float inRad = 0.5f * innerRadius / radius; // Convert inner radius to normalized space

    // Set uniforms
    GLint loc = glGetUniformLocation(shaderProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

    loc = glGetUniformLocation(shaderProgram, "scale");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));

    loc = glGetUniformLocation(shaderProgram, "uColor");
    glUniform4fv(loc, 1, glm::value_ptr(color));

    loc = glGetUniformLocation(shaderProgram, "inRad");
	glUniform1f(loc, inRad); // Set inner radius

    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, locs.size());
    glBindVertexArray(0);


}

void
MeshPointRendering::Cleanup() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (instanceVBO) glDeleteBuffers(1, &instanceVBO);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (shaderProgram) glDeleteProgram(shaderProgram);
}