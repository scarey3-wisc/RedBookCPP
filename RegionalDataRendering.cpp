#include "RegionalDataRendering.h"
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <string>
#include "RenderUtils.h"
#include "ShaderSource.h"

static float rectVertices[] = {
    1.0f, 0.0f,
    0.0f, 0.0f,
    1.0f, 1.0f,
    0.0f, 1.0f
};

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

static const std::string bandedHeightFragSource = 
sourceHeader +
hsb2rgb
+ R"(
in vec2 fragLocalPos; // Rectangle vertex from vertex shader
out vec4 FragColor;

uniform isampler2D uHeightmap;

void main() 
{
    int h = texture(uHeightmap, fragLocalPos).r; // sample red channel
    if(h < 0) h = 0; // clamp negative heights to zero
    if(h == 0)
        discard;
    
    float height = float(h) / 65536.0;
    if(height < 0)
        discard;
    
    float bandElev = mod(500 * pow(height, 1.0 / 3.0), 1000);
	float hue = bandElev / 1000;
	float brt = 0.2 + height / 7000;
    if(brt > 1.0)
        brt = 1.0;
		
	FragColor = vec4(hsb2rgb(vec3(hue, 0.8, brt)), 1.0);
    
}
)";

//------------------------------------------------------------------------------

static const std::string hillshadeFragSource =
sourceHeader +
hsb2rgb +
getHeightmapNormal + R"(

in vec2 fragLocalPos; // Rectangle vertex from vertex shader
out vec4 FragColor;

uniform isampler2D uHeightmap;
uniform float pixelWidthInMeters;
uniform vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));

void main() 
{
    int h = texture(uHeightmap, fragLocalPos).r; // sample red channel
    if(h < 0) h = 0; // clamp negative heights to zero
    if(h == 0)
        discard;
    
    float height = float(h) / 65536.0;
    if(height < 0)
        discard;

    vec3 normal = getHeightmapNormal(uHeightmap, fragLocalPos, pixelWidthInMeters, 65536.0);
    float light = dot(normal, lightDir);
	float brt = (light + 0.1) / 1.1;

    if(brt < 0)
        brt = 0;
    if(brt > 1.0)
        brt = 1.0;
    
    float sat = 0.0;
    float hue = 0.0;
    if(height < 50)
    {
        hue = 56.0 / 360.0;
        sat = .51;
    }
    else if (height < 200)
    {
        hue = 110.0 / 360.0;
        sat = 0.67;
    }
    else if(height < 1000)
    {
        hue = 129.0 / 360.0;
        sat = 0.71;
    }
    else if(height < 3000)
    {
        hue = 38.0 / 360.0;
        sat = 0.27;
    }
    else
    {
        hue = 0.0;
        sat = 0.0;
    }
		
	FragColor = vec4(hsb2rgb(vec3(hue, sat, brt)), 1.0);
    
}
)";

//------------------------------------------------------------------------------

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
    bandedHeightProgram = CompileAndLinkShader(vertexSource.c_str(), bandedHeightFragSource.c_str());
	hillshadeProgram = CompileAndLinkShader(vertexSource.c_str(), hillshadeFragSource.c_str());
}

//------------------------------------------------------------------------------

void
RegionalDataRendering::RenderHillshade(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap, float metersPerPixel)
{

    glUseProgram(hillshadeProgram);


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
    glUniform1i(glGetUniformLocation(hillshadeProgram, "uHeightmap"), 0);

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
    GLint loc = glGetUniformLocation(hillshadeProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

    loc = glGetUniformLocation(hillshadeProgram, "aOffset");
    glUniform2fv(loc, 1, glm::value_ptr(location));

    loc = glGetUniformLocation(hillshadeProgram, "scale");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));

    loc = glGetUniformLocation(hillshadeProgram, "pixelWidthInMeters");
    glUniform1f(loc, metersPerPixel);


    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
    glBindVertexArray(0);

    glDeleteTextures(1, &tex);
}

//------------------------------------------------------------------------------

void
RegionalDataRendering::RenderBandedHeight(glm::vec2 location, float width, float height, const glm::mat4& viewProj, const int32_t* heightmap)
{

    glUseProgram(bandedHeightProgram);


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
    glUniform1i(glGetUniformLocation(bandedHeightProgram, "uHeightmap"), 0);

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
    GLint loc = glGetUniformLocation(bandedHeightProgram, "uViewProj");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mvp));

	loc = glGetUniformLocation(bandedHeightProgram, "aOffset");
	glUniform2fv(loc, 1, glm::value_ptr(location));

    loc = glGetUniformLocation(bandedHeightProgram, "scale");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(scale));


    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
    glBindVertexArray(0);

	glDeleteTextures(1, &tex);
}

//------------------------------------------------------------------------------

void
RegionalDataRendering::Cleanup() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (bandedHeightProgram) glDeleteProgram(bandedHeightProgram);
	if (hillshadeProgram) glDeleteProgram(hillshadeProgram);
}