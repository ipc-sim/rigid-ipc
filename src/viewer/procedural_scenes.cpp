#include "viewer.hpp"

#include <viewer/imgui_ext.hpp>

namespace ccd {

void ViewerMenu::draw_procedural_scene_menu()
{
    static int scene_id = static_cast<int>(ProceduralScene::LINE_STACK);
    ImGui::Combo("scene##procedural", &scene_id, ProceduralSceneNames,
        IM_ARRAYSIZE(ProceduralSceneNames));

    ProceduralScene scene = static_cast<ProceduralScene>(scene_id);
    switch (scene) {
    case ProceduralScene::LINE_STACK:
        draw_line_stack();
        break;
    case ProceduralScene::CHAIN:
        draw_chain_menu();
        break;
    }
}

void ViewerMenu::create_line_stack_scene(
    int num_lines, double displacment_scale)
{
    Eigen::MatrixX2d vertices(2 * num_lines + 2, 2);
    Eigen::MatrixX2d displacements
        = Eigen::MatrixX2d::Zero(2 * num_lines + 2, 2);
    Eigen::MatrixX2i edges(num_lines + 1, 2);

    vertices.row(0) << -0.05, 0.1;
    vertices.row(1) << 0.05, 0.2;
    displacements(0, 1) = -1;
    displacements(1, 1) = -1;
    displacements *= displacment_scale;
    edges.row(0) << 0, 1;

    Eigen::VectorXd ys = Eigen::VectorXd::LinSpaced(num_lines, 0, -1);
    for (int i = 1; i < edges.rows(); i++) {
        vertices.row(2 * i) << -1, ys(i - 1);
        vertices.row(2 * i + 1) << 1, ys(i - 1);
        edges.row(i) << 2 * i, 2 * i + 1;
    }

    state.vertices = vertices;
    state.displacements = displacements;
    state.edges = edges;
    state.reset_scene();
    state_history.push_back(state);
    load_state();
}

void ViewerMenu::draw_line_stack()
{
    static int num_lines = 3;
    static double displacment_scale = 10;
    ImGui::InputIntBounded("line count##line-stack", &num_lines, 0,
        std::numeric_limits<int>::max(), 1, 10);
    ImGui::InputDouble("scale disp.##line-stack", &displacment_scale);
    if (ImGui::Button("Make Line Stack##Edit", ImVec2(-1, 0))) {
        create_line_stack_scene(num_lines, displacment_scale);
    }
}

void ViewerMenu::create_chain_scene(int num_links)
{
    state.load_scene(std::string(FIXTURES_DIR) + "/chain/one-links.json");

    Eigen::MatrixX2d chain_vertices = state.vertices.replicate(num_links, 1);
    Eigen::MatrixX2d chain_displacements
        = state.displacements.replicate(num_links, 1);
    Eigen::MatrixX2i chain_edges = state.edges.replicate(num_links, 1);

    long num_vertices = state.vertices.rows(), num_edges = state.edges.rows();
    for (int i = 0; i < num_links; i++) {
        chain_vertices.block(i * num_vertices, 0, num_vertices, 2)
            .col(1)
            .array()
            -= 2 * i;
        chain_displacements.block(i * num_vertices, 0, num_vertices, 2)
            .col(1)
            .array()
            -= 0.25 * i;
        chain_edges.block(i * num_edges, 0, num_edges, 2).array()
            += i * num_vertices;
    }

    state.vertices = chain_vertices;
    state.displacements = chain_displacements;
    state.edges = chain_edges;

    state.fit_scene_to_canvas();
    state.reset_scene();

    state_history.push_back(state);
    load_state();
}

void ViewerMenu::draw_chain_menu()
{
    static int num_links = 2;
    static double scale_displacment = 10;
    ImGui::InputIntBounded("link count##chain", &num_links, 1,
        std::numeric_limits<int>::max(), 1, 10);
    // ImGui::InputDouble("scale disp.##chain", &scale_displacment);
    if (ImGui::Button("Make Chain##chain", ImVec2(-1, 0))) {
        create_chain_scene(num_links);
    }
}

} // namespace ccd
