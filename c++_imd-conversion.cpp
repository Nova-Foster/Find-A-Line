#include <iostream>
#include <fstream>
#include <vector>

int main() {
    std::ifstream file("test.imd", std::ios::binary);
    if (file.is_open()) {
        unsigned short int VER, W, H;
        file.read(reinterpret_cast<char*>(&VER), sizeof(unsigned short int));
        file.read(reinterpret_cast<char*>(&W), sizeof(unsigned short int));
        file.read(reinterpret_cast<char*>(&H), sizeof(unsigned short int));

        std::vector<signed int> m_s32data(W * H);
        file.read(reinterpret_cast<char*>(m_s32data.data()), sizeof(signed int) * W * H);
        file.close();

        std::vector<float> fdata(W * H);

        for (int k = 0; k < H * W; k++) {
            if (m_s32data[k] != 0) {
                fdata[k] = static_cast<float>(m_s32data[k]) / 1000.0f;
            } else {
                fdata[k] = 0.0f;
            }
        }

        // No need to manually deallocate memory with vectors

    } else {
        std::cerr << "Failed to open the file." << std::endl;
    }
    
    return 0;
}