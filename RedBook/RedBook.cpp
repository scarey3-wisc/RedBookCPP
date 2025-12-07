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
#include "VoronoiAlgorithms.h"
#include "ContinentGenAlgorithms.h"
#include <chrono>
#include <omp.h>
#include "GlobalRand.h"
#include "..\HydrologySolver\FileGlobals.h"

#include "FluidSolver.h"
#include <thread>

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

double renderTimeMovingAverage = 15.0;

constexpr int SCALE_LABEL_TARGET_WIDTH = 200;
constexpr int NUM_LABELS = 18;
constexpr const char* TARGET_SCALE_LABLES[] =
{
    "200 Miles",
    "100 Miles",
    "50 Miles",
    "25 Miles",
    "10 Miles",
    "5 Miles",
    "2 Miles",
    "1 Mile",
    "1/2 Mile",
    "1/3 Mile",
    "1/4 Mile",
    "1000 Feet",
    "500 Feet",
    "200 Feet",
    "100 Feet",
    "50 Feet",
    "20 Feet",
    "10 Feet"
};
constexpr double LABEL_DISTANCES[] =
{
    321869,
    160934,
    80467.2,
    40233.6,
    16093.4,
    8046.72,
    3218.69,
    1609.34,
    804.672,
    536.448,
    406.336,
    304.8,
    152.4,
    60.96,
    30.48,
    15.24,
    6.096,
    3.048
};




WorldMap* myWorldMap = nullptr;


#include "Heightmap.h"

void StartupTests()
{

}

void RedBookInfoPanel(float panelWidth)
{
    ImGui::Text("Render Info");
    ImGui::Separator();
    ImGui::Spacing();
    for (int i = 0; i < REGIONAL_MAP_NUM_LODS; i++)
    {
        int LOD = i + 1;
        int max = myWorldMap->heightmaps.GetMaxCapacity(LOD);
        int stored = myWorldMap->heightmaps.GetNumUsed(LOD);
        int pinned = 0;
		ImGui::Text("LOD #%d: %d / %d (%d)", LOD, stored, max, pinned); // Placeholder for actual cache size
    }
	size_t numQueued = myWorldMap->myWorkerThread.NumActiveRequests();
	ImGui::Text("%d queued", numQueued); // Placeholder for render queue size
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

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	GLuint textureColorbuffer = myWorldMap->Render((int)totalWidth, (int)totalHeight);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    ImVec2 pos = ImGui::GetCursorScreenPos();

    ImGui::Image((void*)(intptr_t)textureColorbuffer, ImVec2(totalWidth, totalHeight),
        ImVec2(0, 1), ImVec2(1, 0));

    double milis = (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    double deltaAverageSmoothness = 0.02;
    renderTimeMovingAverage = (1.0 * (1.0 - deltaAverageSmoothness) * renderTimeMovingAverage + 
        deltaAverageSmoothness * milis);

	int fps = (int) (1000.0 / renderTimeMovingAverage);
	string fpsText = "FPS: " + to_string(fps) + " (" + to_string((int) renderTimeMovingAverage) + " ms)";

    ImDrawList* drawList = ImGui::GetWindowDrawList();


    ImVec2 fpsTextSize = ImGui::CalcTextSize(fpsText.c_str());

    ImVec2 fpsRectEnd(pos.x + fpsTextSize.x + 6, pos.y + fpsTextSize.y + 4);
    drawList->AddRectFilled(pos, fpsRectEnd, IM_COL32(255, 255, 255, 255));

    drawList->AddText(ImVec2(pos.x + 3, pos.y + 2), IM_COL32(0, 0, 0, 255), fpsText.c_str());


    double mPerTile = LocalMap::METER_DIM;
    double mPerPixel = mPerTile / myWorldMap->tileSize;
    double mInTargetScale = SCALE_LABEL_TARGET_WIDTH * mPerPixel;
    double mInChosenScale = -1;
    string chosenScale;
    for (int i = 0; i < NUM_LABELS; i++)
    {
        double mInScale = LABEL_DISTANCES[i];
        string scaleName = TARGET_SCALE_LABLES[i];
        if (i == 0 || (abs(mInTargetScale - mInScale) < abs(mInTargetScale - mInChosenScale)))
        {
            mInChosenScale = mInScale;
            chosenScale = scaleName;
        }
    }
    float barWidth = (float) (mInChosenScale / mPerPixel);

    ImVec2 textSize = ImGui::CalcTextSize(chosenScale.c_str());

    float h = textSize.y + 15;
    float w = max(textSize.x, barWidth) + 10;

    pos.x += totalWidth;
	pos.y += totalHeight;
    pos.x -= w;
    pos.y -= h;
    pos.x -= 10;
    pos.y -= 10;
    pos.x += w / 2;

    ImVec2 rectStart(pos.x - 0.5f * textSize.x - 5, pos.y);
	ImVec2 rectEnd(pos.x + 0.5f * textSize.x + 5, pos.y + textSize.y);
    drawList->AddRectFilled(rectStart, rectEnd, IM_COL32(255, 255, 255, 255));

	ImVec2 textStart(pos.x - 0.5f * textSize.x, pos.y);
    drawList->AddText(textStart, IM_COL32(0, 0, 0, 255), chosenScale.c_str());

    pos.y += h / 2;
    float lineY = h / 2 - 5;

	ImVec2 lineStart(pos.x - 0.5f * barWidth, pos.y + lineY - 2);
	ImVec2 lineEnd(pos.x + 0.5f * barWidth, pos.y + lineY + 2);
    drawList->AddRectFilled(lineStart, lineEnd, IM_COL32(255, 187, 0, 255));

    ImVec2 vertLineStart(lineStart.x - 1, pos.y + lineY - 6);
	ImVec2 vertLineEnd(lineStart.x + 1, pos.y + lineY + 6);
    drawList->AddRectFilled(vertLineStart, vertLineEnd, IM_COL32(255, 187, 0, 255));

    vertLineStart.x = lineEnd.x - 1;
	vertLineEnd.x = lineEnd.x + 1;
    drawList->AddRectFilled(vertLineStart, vertLineEnd, IM_COL32(255, 187, 0, 255));

    vertLineStart.x = pos.x - 1;
	vertLineEnd.x = pos.x + 1;
    vertLineStart.y += 2;
    vertLineEnd.y -= 2;
    drawList->AddRectFilled(vertLineStart, vertLineEnd, IM_COL32(255, 187, 0, 255));
}

void RedBookToolPanel(float panelWidth)
{
    // Fixed panel content
    ImGui::Text("Fixed right panel");
}

void RedBookDisplay(int screenWidth, int screenHeight)
{
    ImGui::SetNextWindowSize(ImVec2((float) screenWidth, (float) screenHeight));
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
    std::array<int, REGIONAL_MAP_NUM_LODS> capacities = std::to_array(REGIONAL_RASTER_CAPCITIES);
    myWorldMap = new WorldMap("Nerean Sea", capacities);
    vector<SamplePoint*> vp;
    myWorldMap->FillAllContinents(0, 0, vp);
    StartupTests();

    vector<SamplePoint*> coastal = VoronoiAlgorithms::FindBoundaryPoints(vp, TerrainTemplate::OCEAN);
    VoronoiAlgorithms::ConvertSeasToLakes(coastal, Switches::MAX_SAMPLE_POINTS_IN_LAKE);
    VoronoiAlgorithms::ConvertCoastalLakeToOcean(vp);

    vector<SamplePoint*> continent = VoronoiAlgorithms::FindAllWithinBoundary(vp, TerrainTemplate::OCEAN);
    ContinentGenAlgorithms::BlurUpliftForTectonicAlgorithm(continent, 25, 5);

    vector<SamplePoint*> rough = VoronoiAlgorithms::FindAllOfType(continent, TerrainTemplate::ROUGH);
    vector<SamplePoint*> mountains = VoronoiAlgorithms::FindAllOfType(continent, TerrainTemplate::MOUNTAINS);

    //paper recommends 2.5 * 10^5 and 5.611 * 10^-7 for erosion
    cout << "Tectonic Uplift Algo: Detail-0" << endl << endl;
    //ContinentGenAlgorithms.RunTectonicUpliftAlgorithm(continent, 2.5 * 100000, 5.611 * 0.0000001, 300, 0.0002);

    for (SamplePoint* sp : rough)
    {
        if (Rand::Float() < 0.99)
            continue;
        sp->MakeLake();
    }

    VoronoiAlgorithms::IncreaseFractureLevel(rough);
    cout << "Tectonic Uplift Algo: Detail-1" << endl << endl;
    //ContinentGenAlgorithms.RunTectonicUpliftAlgorithm(rough, 2.5 * 100000, 5.611 * 0.0000001, 300, 0.0002);

    for (SamplePoint* sp : mountains)
    {
        if (Rand::Float() < 0.96)
            continue;
        sp->MakeLake();
    }

    VoronoiAlgorithms::IncreaseFractureLevel(mountains);
    cout << "Tectonic Uplift Algo: Detail-2" << endl << endl;
    //ContinentGenAlgorithms.RunTectonicUpliftAlgorithm(mountains, 2.5 * 100000, 5.611 * 0.0000001, 300, 0.0002);


    cout << "Tectonic Uplift Algo: River Reset" << endl;
    ContinentGenAlgorithms::RunTectonicUpliftAlgorithm(continent, 2.5 * 100000, 5.611 * 0.0000001, 500, 0.0002);
    cout << endl;

    std::thread t([]() {
        RegionalMap* center = myWorldMap->GetRegion(0, 0);
        RegionalDataLoc top = RegionalDataLoc(center->GetWorldX(), center->GetWorldY(), 0, 0, 1);
        FluidSolver solveSomeFluid(1);
        solveSomeFluid.RecursiveSolutionCycle(top, myWorldMap->heightmaps, myWorldMap->depthmaps, 3);
        //solveSomeFluid.FullSolutionCycle(top, myWorldMap->heightmaps, myWorldMap->depthmaps);
    });
    t.detach();
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

    wchar_t buffer[MAX_PATH];
    DWORD size = GetModuleFileNameW(NULL, buffer, MAX_PATH);
    std::filesystem::path exePath = std::filesystem::path(buffer, buffer + size);
    std::filesystem::path exeDir = exePath.parent_path();
    std::filesystem::path topDir = exeDir.parent_path().parent_path().parent_path();
    std::filesystem::path outPath = topDir;
    MY_PATH = outPath;
    MY_PATH.append("output");

	omp_set_num_threads(omp_get_max_threads() * 2 / 3);
    AllocConsole();
    FILE* fp;
    freopen_s(&fp, "CONOUT$", "w", stdout);
    freopen_s(&fp, "CONOUT$", "w", stderr);

    std::cout << "Using " << omp_get_max_threads() << " threads for parallel processing." << std::endl;


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
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
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
    myWorldMap->myWorkerThread.DiscardAllWaitingTasks();
    myWorldMap->myWorkerThread.TerminateThreads();

    delete myWorldMap;
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
