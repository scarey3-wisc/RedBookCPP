#pragma once
#include <random>
class Rand {
public:

    static int Int(int min, int max) 
    {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(engine());
    }

    static float Float(float min = 0.0f, float max = 1.0f) 
    {
        std::uniform_real_distribution<float> dist(min, max);
        return dist(engine());
    }

private:
    static std::mt19937& engine()
    {
        static std::mt19937 rng{ std::random_device{}() }; // seeded once
        return rng;
    }
};