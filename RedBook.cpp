#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#include <glad/glad.h>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h> // Will drag system OpenGL headers
#include <Windows.h>
#include "RedBook.h"
#include "RedBookGlobal.h"
#include "Switches.h"
#include <iostream>
#include "WorldMap.h"
#include "SamplePoint.h"
#include <algorithm>

using namespace std;
using namespace RedBook;

ImVec4 g_WhiteColor = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
ImVec4 g_BlackColor = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
ImVec4 g_GrayBackgroundColor = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);
ImVec4 g_LightBlueTextColor = ImVec4(0.8f, 0.87f, 1.0f, 1.0f);
ImVec4 g_LightBrownTextColor = ImVec4(0.8f, 0.7f, 0.6f, 1.0f);

ImVec4 g_LightGreyButtonColor = ImVec4(0.4f, 0.4f, 0.4f, 1.0f);
ImVec4 g_MediumGreyButtonColor = ImVec4(0.3f, 0.3f, 0.3f, 1.0f);
ImVec4 g_DarkGreyButtonColor = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);

WorldMap* myWorldMap = nullptr;


void RedBookInfoPanel(float panelWidth)
{
    ImGui::Text("Render Info");
    ImGui::Separator();
    ImGui::Spacing();
    for (int i = 0; i < TILE_RENDERING_CACHES; i++)
    {
		ImGui::Text("%s: %ld", TILE_RENDERING_NAMES[i], 50); // Placeholder for actual cache size
    }
	ImGui::Text("%d queued (placeholder)", 10); // Placeholder for render queue size
    ImGui::NewLine();

	int dimRange[] = { 0, 1, 2, 3, 4, 5 }; // Placeholder for actual dimension range
    for (int i = 0; i < sizeof(dimRange) / sizeof(dimRange[0]); i++)
    {
        ImGui::Text("Dimension %d: %d", i, dimRange[i]); // Placeholder for actual dimension data
	}
	ImGui::Text("Smart Image Count: %d", 100); // Placeholder for smart image count
    ImGui::NewLine();
    ImGui::NewLine();


    ImGui::PushStyleColor(ImGuiCol_Button, g_LightGreyButtonColor); // normal
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, g_MediumGreyButtonColor); // hover
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, g_DarkGreyButtonColor); // pressed
    ImGui::PushStyleColor(ImGuiCol_FrameBg, g_LightGreyButtonColor);      // box background
    ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, g_MediumGreyButtonColor); // hovered background
    ImGui::PushStyleColor(ImGuiCol_FrameBgActive, g_DarkGreyButtonColor);  // active background
    ImGui::PushStyleColor(ImGuiCol_CheckMark, g_LightBrownTextColor);

    ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 6.0f);

	ImGui::Text("Render Options");
    ImGui::Separator();
    ImGui::Spacing();
    if (ImGui::Button("Clear Image Cache"))
    {
        // Clear image cache logic here
    }
	ImGui::Checkbox("Toggle Outlines", &Switches::OUTLINE_MAPS);
	ImGui::Checkbox("Toggle Voronoi Centers", &Switches::PAINT_VORONOI_CENTERS);

	ImGui::NewLine();
    ImGui::NewLine();
    ImGui::Text("Display Paint Modes");
    ImGui::Separator();
    ImGui::Spacing();
    if (ImGui::RadioButton("Paint Terrain", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::TERRAIN_PRETTY))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::TERRAIN_PRETTY;
    if (ImGui::RadioButton("Paint Photos", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::PHOTOGRAPHY))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::PHOTOGRAPHY;

    ImGui::NewLine();
	ImGui::NewLine();
    ImGui::Text("Informational Paint Modes");
    ImGui::Separator();
	ImGui::Spacing();
    if (ImGui::RadioButton("Paint Heightmap", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::ELEVATION_CURR))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::ELEVATION_CURR;
    if (ImGui::RadioButton("Paint Terrain Type", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::TERRAIN))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::TERRAIN;
    if (ImGui::RadioButton("Paint Contour Map", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::CONTOUR))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::CONTOUR;
    if (ImGui::RadioButton("Paint Topography", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::ELEV_GRADIENT))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::ELEV_GRADIENT;

    ImGui::NewLine();
    ImGui::NewLine();
    ImGui::Text("Diagnostic Paint Modes");
    ImGui::Separator();
    ImGui::Spacing();
    if (ImGui::RadioButton("Paint Perlin Noise", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::PERLIN_DISPLAY))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::PERLIN_DISPLAY;
    if (ImGui::RadioButton("Paint Voronoi Cells", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::VORONOI_PURE))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::VORONOI_PURE;
    if (ImGui::RadioButton("Paint Voronoi Interpolation", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::VORONOI_INTERPOLATED))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::VORONOI_INTERPOLATED;
    if (ImGui::RadioButton("Paint Voronoi Triangles", Switches::CURR_PAINT_TYPE == Switches::PAINT_TYPE::VORONOI_TRIANGLES))
        Switches::CURR_PAINT_TYPE = Switches::PAINT_TYPE::VORONOI_TRIANGLES;

	ImGui::PopStyleColor(7); // Pop the button styles
	ImGui::PopStyleVar(1); // Pop the frame rounding style
}

void RedBookDisplayPanel()
{


    float totalWidth = ImGui::GetContentRegionAvail().x;
    float totalHeight = ImGui::GetContentRegionAvail().y;
    // Detect if the mouse is over this specific child
    if (ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem)) {
        // Only pan if the left mouse button is held down
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
            ImVec2 delta = ImGui::GetIO().MouseDelta;

            // Apply to your camera/view offsets
            myWorldMap->dX += delta.x;
            myWorldMap->dY -= delta.y;
        }

        float scroll = ImGui::GetIO().MouseWheel;
        if (scroll != 0.0f) {
            ImVec2 mouse = ImGui::GetMousePos();
            ImVec2 windowPos = ImGui::GetWindowPos();
            ImVec2 windowCursor(mouse.x - windowPos.x, mouse.y - windowPos.y);

			int oldTileSize = myWorldMap->tileSize;

            float zoomFactor = (scroll > 0 ? 1.1f : 0.9f);
            int newTileSize = clamp(int(oldTileSize * zoomFactor + 0.5), WorldMap::MINIMUM_TILE_SIZE, WorldMap::MAXIMUM_TILE_SIZE);

            float ratio = float(newTileSize) / float(oldTileSize);

            windowCursor.x -= totalWidth / 2;
            windowCursor.y -= totalHeight / 2;
            windowCursor.y *= -1;
            double newDX = myWorldMap->dX * ratio - windowCursor.x * (ratio - 1);
            double newDY = myWorldMap->dY * ratio - windowCursor.y * (ratio - 1);
            myWorldMap->dX = newDX;
            myWorldMap->dY = newDY;

            myWorldMap->tileSize = newTileSize;
        }
    }

    if( myWorldMap == nullptr)
    {
        ImGui::Text("World map not initialized.");
        return;
	}
	GLuint textureColorbuffer = myWorldMap->Render((int)totalWidth, (int)totalHeight);
    ImGui::Image((void*)(intptr_t)textureColorbuffer, ImVec2(totalWidth, totalHeight),
        ImVec2(0, 1), ImVec2(1, 0));
}

void RedBookToolPanel(float panelWidth)
{
    // Fixed panel content
    ImGui::Text("Fixed right panel");
}

void RedBookDisplay(int screenWidth, int screenHeight)
{
    ImGui::SetNextWindowSize(ImVec2(screenWidth, screenHeight));
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::Begin("Main Panel", nullptr, 
        ImGuiWindowFlags_NoCollapse | 
        ImGuiWindowFlags_NoMove | 
        ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoDecoration
    );
    float totalWidth = ImGui::GetContentRegionAvail().x;
    float totalHeight = ImGui::GetContentRegionAvail().y;

    ImGui::PushStyleColor(ImGuiCol_Text, g_LightBrownTextColor);
    ImGui::PushStyleColor(ImGuiCol_ChildBg, g_GrayBackgroundColor);
    // ----- Panel 1: Scrollable, auto-width -----
    ImGui::BeginChild("InfoPanel", ImVec2(300, totalHeight), true);

    RedBookInfoPanel(300.0f);
    // Save computed width
    float leftPanelWidth = ImGui::GetWindowSize().x;
    ImGui::EndChild();
	ImGui::PopStyleColor(2); // Pop the styles we pushed

    // Line up next panel horizontally
    ImGui::SameLine();

    float rightPanelWidth = 300.0f;
    float spacing = ImGui::GetStyle().ItemSpacing.x;
    float middlePanelWidth = totalWidth - leftPanelWidth - rightPanelWidth - spacing * 2;

    ImGui::BeginChild("DisplayPanel", ImVec2(middlePanelWidth, totalHeight), ImGuiChildFlags_None);
	RedBookDisplayPanel();
    ImGui::EndChild();

    ImGui::SameLine();

    ImGui::PushStyleColor(ImGuiCol_Text, g_LightBrownTextColor);
    ImGui::PushStyleColor(ImGuiCol_ChildBg, g_GrayBackgroundColor);
    ImGui::BeginChild("ToolPanel", ImVec2(rightPanelWidth, totalHeight), true);
	RedBookToolPanel(rightPanelWidth);
    ImGui::EndChild();
	ImGui::PopStyleColor(2); // Pop the styles we pushed

    ImGui::End();
}

void 
RedBookInitWorld()
{
    myWorldMap = new WorldMap("Nerean Sea");
    vector<SamplePoint*> vp;
    myWorldMap->FillAllContinents(0, 0, vp);
}


static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

// Main code
int WINAPI WinMain(
    HINSTANCE hInstance, 
    HINSTANCE hPrevInstance,
    LPSTR lpCmdLine, 
    int nCmdShow)
{
    AllocConsole();
    FILE* fp;
    freopen_s(&fp, "CONOUT$", "w", stdout);
    freopen_s(&fp, "CONOUT$", "w", stderr);

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
    // GL ES 2.0 + GLSL 100 (WebGL 1.0)
    const char* glsl_version = "#version 100";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(IMGUI_IMPL_OPENGL_ES3)
    // GL ES 3.0 + GLSL 300 es (WebGL 2.0)
    const char* glsl_version = "#version 300 es";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
    // GL 3.2 + GLSL 150
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
    // GL 3.0 + GLSL 130
    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

    // Create window with graphics context
    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
    float main_scale = ImGui_ImplGlfw_GetContentScaleForMonitor(glfwGetPrimaryMonitor()); // Valid on GLFW 3.3+ only
    GLFWwindow* window = glfwCreateWindow((int)(1280 * main_scale), (int)(800 * main_scale), "Red Book World Manager", nullptr, nullptr);
    if (window == nullptr)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync


    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD\n";
        //return -1;
    }

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup scaling
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);        // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
    style.FontScaleDpi = main_scale;        // Set initial font scale. (using io.ConfigDpiScaleFonts=true makes this unnecessary. We leave both here for documentation purpose)

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifdef __EMSCRIPTEN__
    ImGui_ImplGlfw_InstallEmscriptenCallbacks(window, "#canvas");
#endif
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    // - Our Emscripten build process allows embedding fonts to be accessible at runtime from the "fonts/" folder. See Makefile.emscripten for details.
    //style.FontSizeBase = 20.0f;
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf");
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf");
    //IM_ASSERT(font != nullptr);

    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    RedBookInitWorld();
    while (!glfwWindowShouldClose(window))
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        glfwPollEvents();
        if (glfwGetWindowAttrib(window, GLFW_ICONIFIED) != 0)
        {
            ImGui_ImplGlfw_Sleep(10);
            continue;
        }
        
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);

        RedBookDisplay(fbWidth, fbHeight);
        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    delete myWorldMap;
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
