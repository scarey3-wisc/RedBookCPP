#include "RegionalDataRendering.h"
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <string>
#include "RenderUtils.h"

static float rectVertices[] = {
    1.0f, 0.0f,
    0.0f, 0.0f,
    1.0f, 1.0f,
    0.0f, 1.0f
};

static const std::string vertexSource = R"(
#version 330 core

layout (location = 0) in vec2 aPos;      // Rectangle vertex

uniform vec2 aOffset;   // Instance position

uniform mat4 uViewProj; // Your orthographic projection
uniform mat4 scale;

out vec2 fragLocalPos;  // Pass to fragment shader

void main()
{
    gl_Position = uViewProj * (scale * vec4(aPos, 0.0, 1.0) + vec4(aOffset, 0.0, 0.0));
    fragLocalPos = aPos;
}
)";

static const std::string heightmapFragSource = R"(
#version 330 core

in vec2 fragLocalPos; // Rectangle vertex from vertex shader
out vec4 FragColor;

uniform isampler2D uHeightmap;
uniform vec4 uColor;

vec3 hsb2rgb(in vec3 c) {
    vec3 rgb = clamp(
        abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0,
        0.0,
        1.0
    );
    rgb = rgb * rgb * (3.0 - 2.0 * rgb); // smooth interpolation
    return c.z * mix(vec3(1.0), rgb, c.y);
}

void main() 
{
    int h = texture(uHeightmap, fragLocalPos).r; // sample red channel
    if(h < 0) h = 0; // clamp negative heights to zero
    if(h == 0)
        discard;
    
    float height = float(h) / 65535.0;
    if(height < 0)
        discard;
    
    float bandElev = mod(500 * pow(height, 1.0 / 3.0), 1000);
	float hue = bandElev / 1000;
	float brt = 0.2 + height / 10000;
		
	FragColor = vec4(hsb2rgb(vec3(hue, 0.8, brt)), 1.0);
    
}
)";

void
RegionalDataRendering::Init() {
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

    glBindVertexArray(0);

    // Compile shaders (vertex & fragment)
    // You'll need to write helper functions to compile/link shaders
    shaderProgram = CompileAndLinkShader(vertexSource.c_str(), heightmapFragSource.c_str());
}

void
RegionalDataRendering::Render(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap, const glm::vec4& color)
{

    glUseProgram(shaderProgram);


    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    // Upload as 32-bit signed integer texture
    glTexImage2D(
        GL_TEXTURE_2D,
        0,
        GL_R32I,                // internal format (signed 32-bit ints)
        256,
        256,
        0,
        GL_RED_INTEGER,         // format: red channel only, integer interpretation
        GL_INT,                 // type: 32-bit signed integer
        heightmap
    );

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(glGetUniformLocation(shaderProgram, "uHeightmap"), 0);

    // Wrapping control
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Filtering (must be nearest for integer textures)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // Build model matrix: translate and scale unit rectangle
    glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(width, height, 1.0f));

    glm::mat4 mvp = viewProj;

    // Set uniforms
    GLint loc = glGetUniformLocation(shaderProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

	loc = glGetUniformLocation(shaderProgram, "aOffset");
	glUniform2fv(loc, 1, glm::value_ptr(location));

    loc = glGetUniformLocation(shaderProgram, "scale");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));

    loc = glGetUniformLocation(shaderProgram, "uColor");
    glUniform4fv(loc, 1, glm::value_ptr(color));

    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
    glBindVertexArray(0);

	glDeleteTextures(1, &tex);
}

void
RegionalDataRendering::Cleanup() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (shaderProgram) glDeleteProgram(shaderProgram);
}