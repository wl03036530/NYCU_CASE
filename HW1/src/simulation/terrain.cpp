#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Sphere:
            return std::make_unique<SphereTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        case simulation::TerrainType::TiltedPlane:
            return std::make_unique<TiltedPlaneTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    
    // TODO
    Eigen::Vector3f N = this->normal / this->normal.norm();
    for (int i = 0; i < cube.getParticleNum(); i += 1)
    {
        if (abs(N.dot(cube.getParticle(i).getPosition() - this->position)) < eEPSILON) // on the wall
        {
            if (N.dot(cube.getParticle(i).getVelocity()) < 0)  // collision
            {
                Eigen::Vector3f Vn = cube.getParticle(i).getVelocity().dot(N) * N;
                Eigen::Vector3f Vt = cube.getParticle(i).getVelocity() - Vn;
                Eigen::Vector3f V_after = -1 * coefResist * Vn + Vt;
                cube.getParticle(i).setVelocity(V_after);

                Eigen::Vector3f Fc = -(N.dot(cube.getParticle(i).getForce()) * N);
                cube.getParticle(i).addForce(Fc);

                Eigen::Vector3f Ff = -coefFriction * (-N.dot(cube.getParticle(i).getForce()) * Vt);
                cube.getParticle(i).addForce(Ff);
            }
        }
    }
}

// SphereTerrain //

SphereTerrain::SphereTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType SphereTerrain::getType() { return TerrainType::Sphere; }

void SphereTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO
    for (int i = 0; i < cube.getParticleNum(); i += 1)
    {
        Particle p = cube.getParticle(i);
        if ((p.getPosition() - this->position).norm() - this->radius < eEPSILON) // collision
        {
            Eigen::Vector3f N = (cube.getParticle(i).getPosition() - this->position) / (cube.getParticle(i).getPosition() - this->position).norm();
            if (N.dot(cube.getParticle(i).getVelocity()) < 0)  // collision
            {
                Eigen::Vector3f Vn = cube.getParticle(i).getVelocity().dot(N) * N;
                Eigen::Vector3f Vt = cube.getParticle(i).getVelocity() - Vn;
                Eigen::Vector3f V_after = Vn * (cube.getParticle(i).getMass() - this->mass) / (cube.getParticle(i).getMass() + this->mass) + Vt;
                cube.getParticle(i).setVelocity(V_after);

                Eigen::Vector3f Fc = -(N.dot(cube.getParticle(i).getForce()) * N);
                cube.getParticle(i).addForce(Fc);

                Eigen::Vector3f Ff = -coefFriction * (-N.dot(cube.getParticle(i).getForce()) * Vt);
                cube.getParticle(i).addForce(Ff);
            }
        }
    }
}

// BowlTerrain //

BowlTerrain::BowlTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO
    for (int i = 0; i < cube.getParticleNum(); i += 1)
    {
        Particle p = cube.getParticle(i);
        if (this->radius - (p.getPosition() - this->position).norm() <= eEPSILON) // collision
        {
            Eigen::Vector3f N = (cube.getParticle(i).getPosition() - this->position) / (cube.getParticle(i).getPosition() - this->position).norm();
            if (N.dot(cube.getParticle(i).getVelocity()) > 0)  // collision
            {
                Eigen::Vector3f Vn = cube.getParticle(i).getVelocity().dot(N) * N;
                Eigen::Vector3f Vt = cube.getParticle(i).getVelocity() - Vn;
                Eigen::Vector3f V_after = Vn * (cube.getParticle(i).getMass() - this->mass) / (cube.getParticle(i).getMass() + this->mass) + Vt;
                cube.getParticle(i).setVelocity(V_after);

                Eigen::Vector3f Fc = -(N.dot(cube.getParticle(i).getForce()) * N);
                cube.getParticle(i).addForce(Fc);

                Eigen::Vector3f Ff = -coefFriction * (-N.dot(cube.getParticle(i).getForce()) * Vt);
                cube.getParticle(i).addForce(Ff);
            }
        }
    }
}

// TiltedPlaneTerrain //

TiltedPlaneTerrain::TiltedPlaneTerrain() { modelMatrix = util::rotateDegree(0, 0, -45) * util::scale(60, 1, 60); }

TerrainType TiltedPlaneTerrain::getType() { return TerrainType::TiltedPlane; }

void TiltedPlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO
    Eigen::Vector3f N = this->normal / this->normal.norm();
    for (int i = 0; i < cube.getParticleNum(); i += 1)
    {
        if (abs(N.dot(cube.getParticle(i).getPosition() - this->position)) < eEPSILON) // on the wall
        {
            if (N.dot(cube.getParticle(i).getVelocity()) < 0) //collision
            {
                Eigen::Vector3f Vn = cube.getParticle(i).getVelocity().dot(N) * N;
                Eigen::Vector3f Vt = cube.getParticle(i).getVelocity() - Vn;
                Eigen::Vector3f V_after = -1 * coefResist * Vn + Vt;
                cube.getParticle(i).setVelocity(V_after);

                Eigen::Vector3f Fc = -(N.dot(cube.getParticle(i).getForce()) * N);
                cube.getParticle(i).addForce(Fc);

                Eigen::Vector3f Ff = -coefFriction * (-N.dot(cube.getParticle(i).getForce()) * Vt);
                cube.getParticle(i).addForce(Ff);
            }
        }
    }
}
}  // namespace simulation
