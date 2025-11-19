#include "OutlineRendering.h"
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include "RenderUtils.h"

static float rectVertices[] = {
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    0.0f, 1.0f
};

static const char* vertexSource = R"(
#version 330 core

layout (location = 0) in vec2 aPos;      // Rectangle vertex
layout (location = 1) in vec2 aOffset;   // Instance position

uniform mat4 uViewProj; // Your orthographic projection
uniform mat4 scale;

void main()
{
    gl_Position = uViewProj * (scale * vec4(aPos, 0.0, 1.0) + vec4(aOffset, 0.0, 0.0));
}
)";

static const char* fragmentSource = R"(
#version 330 core
out vec4 FragColor;

uniform vec4 uColor;

void main() {
    FragColor = uColor;
}
)";

//------------------------------------------------------------------------------

void 
OutlineRendering::Init() {
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

//------------------------------------------------------------------------------

void 
OutlineRendering::Render(std::vector<glm::vec2>& locs, float width, float height, const glm::mat4& viewProj, const glm::vec4& color)
{

    glUseProgram(shaderProgram);

    // Build model matrix: translate and scale unit rectangle
	glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(width, height, 1.0f));
    

    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, locs.size() * sizeof(glm::vec2), locs.data(), GL_DYNAMIC_DRAW);

    glm::mat4 mvp = viewProj;

    // Set uniforms
    GLint loc = glGetUniformLocation(shaderProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

    loc = glGetUniformLocation(shaderProgram, "scale");
	glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));

    loc = glGetUniformLocation(shaderProgram, "uColor");
    glUniform4fv(loc, 1, glm::value_ptr(color));

    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_LINE_LOOP, 0, 4, locs.size());
    glBindVertexArray(0);
    
   
}

//------------------------------------------------------------------------------

void 
OutlineRendering::Cleanup() {
    if (vbo) glDeleteBuffers(1, &vbo);
	if (instanceVBO) glDeleteBuffers(1, &instanceVBO);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (shaderProgram) glDeleteProgram(shaderProgram);
}