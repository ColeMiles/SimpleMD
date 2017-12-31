#include "simpleMD.hpp"

int main() {
    vec3 box_dim = {1.0, 1.0, 1.0};
    SimpleMDBox box(box_dim, 100, 1.0, 0.0005);
    return 0;
}
